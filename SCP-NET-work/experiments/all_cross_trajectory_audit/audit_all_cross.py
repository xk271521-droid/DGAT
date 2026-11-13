import sys, json, os
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import anndata as ad
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

BASE = Path(r"E:/DGAT")
DATA_ROOT = BASE / "ST2P/data/processed/DGAT_final10_common9384_muon017_axis0_20261101"
RUNS_BASE = BASE / "SCP-NET-work/runs/all_cross_trajectory_runs"
OUTPUT_DIR = BASE / "SCP-NET-work/experiments/all_cross_trajectory_audit"

sys.path.insert(0, str(BASE / "ST2P/src"))
from triscale_center_protein_v9_dgat_compatible import dgat_metrics
from v9_robust_source_subset_test_oracle import corrected_knn, dense, panel
from triscale_center_protein_v2_relpos import relative_position

import importlib.util
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

manifest = json.loads((DATA_ROOT / "prepared_manifest.json").read_text())
names = [entry["dataset"] for entry in manifest["datasets"]]
genes, proteins = panel(DATA_ROOT, names, "rna_dgat.h5ad"), panel(DATA_ROOT, names, "adt_dgat_clr.h5ad")

cognate_map_dict = {
    'CD163': 'CD163', 'CR2': 'CR2', 'PCNA': 'PCNA', 'VIM': 'VIM', 'KRT5': 'KRT5',
    'CD68': 'CD68', 'CEACAM8': 'CEACAM8', 'PTPRC_1': 'PTPRC', 'PTPRC_2': 'PTPRC',
    'PAX5': 'PAX5', 'SDC1': 'SDC1', 'CD8A': 'CD8A', 'BCL2': 'BCL2', 'CD19': 'CD19',
    'PDCD1': 'PDCD1', 'ACTA2': 'ACTA2', 'FCGR3A': 'FCGR3A', 'ITGAX': 'ITGAX',
    'CXCR5': 'CXCR5', 'EPCAM': 'EPCAM', 'MS4A1': 'MS4A1', 'CD3E': 'CD3E',
    'CD14': 'CD14', 'CD40': 'CD40', 'PECAM1': 'PECAM1', 'CD4': 'CD4',
    'ITGAM': 'ITGAM', 'CD27': 'CD27', 'CCR7': 'CCR7', 'CD274': 'CD274',
}
gene_to_idx = {g: i for i, g in enumerate(genes)}
cognate_indices = [gene_to_idx[cognate_map_dict[p]] if p in cognate_map_dict and cognate_map_dict[p] in gene_to_idx else 0 for p in proteins]
has_cognate = [1.0 if p in cognate_map_dict and cognate_map_dict[p] in gene_to_idx else 0.0 for p in proteins]

# List of 10 cross targets
targets_info = [
    ('cross_holdout_LN_A1_LN', 'A1_LN', RUNS_BASE / 'cross_holdout_LN_A1_LN'),
    ('cross_holdout_LN_D1_LN', 'D1_LN', RUNS_BASE / 'cross_holdout_LN_D1_LN'),
    ('cross_holdout_Tonsil_A1_TNSL', 'A1_TNSL', RUNS_BASE / 'cross_holdout_Tonsil_A1_TNSL'),
    ('cross_holdout_Tonsil_D1_TNSL', 'D1_TNSL', RUNS_BASE / 'cross_holdout_Tonsil_D1_TNSL'),
    ('cross_holdout_Tonsil_HumanTonsil_10x', 'HumanTonsil_10x', RUNS_BASE / 'cross_holdout_Tonsil_HumanTonsil_10x'),
    ('cross_holdout_Tonsil_HumanTonsil_AddOns_S1', 'HumanTonsil_AddOns_S1', RUNS_BASE / 'cross_holdout_Tonsil_HumanTonsil_AddOns_S1'),
    ('cross_holdout_Breast_Cancer_Breast_Cancer', 'Breast_Cancer', RUNS_BASE / 'cross_holdout_Breast_Cancer_Breast_Cancer'),
    ('cross_holdout_Breast_Cancer_PBC_PR_6837', 'PBC_PR_6837', RUNS_BASE / 'cross_holdout_Breast_Cancer_PBC_PR_6837'),
    ('cross_holdout_Breast_Cancer_PBC_PR_6835_5A', 'PBC_PR_6835_5A', RUNS_BASE / 'cross_holdout_Breast_Cancer_PBC_PR_6835_5A'),
    ('cross_holdout_GBM', 'GBM', BASE / 'SCP-NET-work/runs/epoch_trajectory_probe/cross_holdout_GBM_all_epochs')
]

summary_rows = []

for d_name, t_name, run_dir in targets_info:
    check_dir = run_dir / 'epoch_checkpoints'
    if not check_dir.exists():
        print(f"Skipping {d_name}, checkpoints not ready yet...")
        continue
    
    print(f">>> Auditing 20 epochs for {t_name} ({d_name})...")
    spec = importlib.util.spec_from_file_location(f"model_{t_name}", run_dir / "script_snapshot.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    config = json.loads((run_dir / "config.json").read_text())
    sources = config["source_sections"]
    
    # Target data
    rna = ad.read_h5ad(DATA_ROOT / t_name / "rna_dgat.h5ad")[:, genes].copy()
    coords = np.asarray(rna.obsm["spatial"], dtype=np.float32)
    neighbours = corrected_knn(coords, config["regional_neighbors"])
    local, regional = neighbours[:, :config["local_neighbors"]], neighbours[:, config["local_neighbors"]:]
    nearest = np.linalg.norm(coords[local[:, 0]] - coords, axis=1)
    scale = float(np.median(nearest[nearest > 0]))
    lp = relative_position(coords, local, scale)
    rp = relative_position(coords, regional, scale)
    x_mat = dense(rna.X)
    adt = ad.read_h5ad(DATA_ROOT / t_name / "adt_dgat_clr.h5ad")[:, proteins]
    target_truth = dense(adt.X)
    target_idx = np.arange(len(x_mat), dtype=np.int64)
    target_section_ids = np.full(len(target_idx), len(sources), dtype=np.int64)

    # Context
    ctx_npz = np.load(run_dir / "rna_slice_context_descriptors.npz")
    canon_npz = np.load(run_dir / "rna_slice_canonicalization.npz")
    desc_t = torch.from_numpy(ctx_npz["descriptors"]).to(device)
    means_t = torch.from_numpy(canon_npz["mean"]).to(device)
    stds_t = torch.from_numpy(canon_npz["std"]).to(device)

    # Calibration parameters
    with open(BASE / f"SCP-NET-work/runs/protein_convex_rna_calibration_10_dataset_20261103/all_directions/{d_name}/calibration_results.json", "r") as f:
        calib_cfg = json.load(f)
    rho = calib_cfg["selected_rho"]
    lambdas = np.asarray([calib_cfg["selected_lambda_by_protein"][p] for p in proteins])

    # Target graph weights
    distances, n_nbrs = NearestNeighbors(n_neighbors=9).fit(coords).kneighbors(coords)
    distances, n_nbrs = distances[:, 1:], n_nbrs[:, 1:]
    n_comp = min(32, x_mat.shape[1], len(x_mat) - 1)
    z_pca = PCA(n_components=n_comp, random_state=42).fit_transform(x_mat)
    f_dist = np.linalg.norm(z_pca[:, None, :] - z_pca[n_nbrs], axis=-1)
    sp_scale = np.median(distances[distances > 0])
    ft_scale = np.median(f_dist[f_dist > 0])
    sp_w = np.exp(-0.5 * (distances / (sp_scale + 1e-8)) ** 2)
    ft_w = np.exp(-0.5 * (f_dist / (ft_scale + 1e-8)) ** 2)
    sp_w /= sp_w.sum(axis=1, keepdims=True)
    ft_w /= ft_w.sum(axis=1, keepdims=True)
    target_weights = rho * sp_w + (1.0 - rho) * ft_w

    model = module.SliceConditionedNicheProtein(
        len(genes), len(proteins), config["hidden_dim"], config["heads"], config["ffn_dim"],
        config["dropout"], config["protein_embedding_dim"], config["router_rank"],
        context_input_dim=desc_t.shape[1], context_rank=config["context_rank"],
        protein_specific_router=not config["disable_protein_specific_router"],
        protein_aware_attention=not config["disable_protein_aware_attention"],
        use_slice_context=not config["disable_slice_context"],
        use_canonicalization=not config["disable_canonicalization"],
        use_context_adapter=not config["disable_context_adapter"],
        cognate_indices=cognate_indices, has_cognate=has_cognate,
        use_cognate_anchor=not config["disable_cognate_anchor"],
    ).to(device)

    hist_df = pd.read_csv(run_dir / "metrics_history.csv")
    orig_best_ep = int(json.loads((run_dir / "results.json").read_text())["best_epoch"])

    records = []
    for ep in range(1, 21):
        pth = check_dir / f"epoch_{ep:02d}.pth"
        model.load_state_dict(torch.load(pth, map_location=device, weights_only=True))
        model.eval()
        with torch.no_grad():
            preds = []
            for start in range(0, len(target_idx), config["batch_size"]):
                idx = target_idx[start:start + config["batch_size"]]
                inp = module.model_batch_inputs(x_mat, local, regional, lp, rp, idx, target_section_ids, desc_t, means_t, stds_t, device)
                preds.append(model(*inp, need_attention=False)["final"].cpu().numpy())
        p_raw = np.concatenate(preds)
        p_cal = (1.0 - lambdas[None, :]) * p_raw + lambdas[None, :] * (target_weights[:, :, None] * p_raw[n_nbrs]).sum(1)
        
        m_cal = dgat_metrics(p_cal, target_truth, proteins)
        records.append({
            "epoch": ep,
            "calib_PCC": m_cal["pcc_pro_mean"],
            "calib_RMSE": m_cal["rmse_global"],
            "calib_MAE": m_cal["mae_global"]
        })
    
    rec_df = pd.DataFrame(records)
    rec_df.to_csv(OUTPUT_DIR / f"{d_name}_trajectory.csv", index=False)

    best_pcc_row = rec_df.loc[rec_df["calib_PCC"].idxmax()]
    best_rmse_row = rec_df.loc[rec_df["calib_RMSE"].idxmin()]
    orig_row = rec_df.loc[rec_df["epoch"] == orig_best_ep].iloc[0]

    summary_rows.append({
        "direction": d_name,
        "target": t_name,
        "orig_epoch": orig_best_ep,
        "orig_PCC": orig_row["calib_PCC"],
        "orig_RMSE": orig_row["calib_RMSE"],
        "best_pcc_epoch": int(best_pcc_row["epoch"]),
        "best_PCC": best_pcc_row["calib_PCC"],
        "best_pcc_epoch_RMSE": best_pcc_row["calib_RMSE"],
        "best_rmse_epoch": int(best_rmse_row["epoch"]),
        "best_RMSE": best_rmse_row["calib_RMSE"],
        "best_rmse_epoch_PCC": best_rmse_row["calib_PCC"]
    })

sum_df = pd.DataFrame(summary_rows)
sum_df.to_csv(OUTPUT_DIR / "CROSS_HOLDOUT_ORACLE_AUDIT_SUMMARY.csv", index=False)
print(">>> Audit finished and summary saved!")
