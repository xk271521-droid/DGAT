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
RUN_DIR = BASE / "SCP-NET-work/runs/epoch_trajectory_probe/cross_holdout_GBM_all_epochs"
CHECKPOINT_DIR = RUN_DIR / "epoch_checkpoints"

sys.path.insert(0, str(BASE / "ST2P/src"))
from triscale_center_protein_v9_dgat_compatible import dgat_metrics
from v9_robust_source_subset_test_oracle import corrected_knn, dense, panel
from triscale_center_protein_v2_relpos import relative_position

# Load model snapshot
import importlib.util
spec = importlib.util.spec_from_file_location("probe_model", RUN_DIR / "script_snapshot.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Load config and panel
config = json.loads((RUN_DIR / "config.json").read_text())
manifest = json.loads((DATA_ROOT / "prepared_manifest.json").read_text())
names = [entry["dataset"] for entry in manifest["datasets"]]
genes, proteins = panel(DATA_ROOT, names, "rna_dgat.h5ad"), panel(DATA_ROOT, names, "adt_dgat_clr.h5ad")
sources = config["source_sections"]

# Load target GBM data
rna = ad.read_h5ad(DATA_ROOT / "GBM/rna_dgat.h5ad")[:, genes].copy()
coords = np.asarray(rna.obsm["spatial"], dtype=np.float32)
neighbours = corrected_knn(coords, config["regional_neighbors"])
local, regional = neighbours[:, :config["local_neighbors"]], neighbours[:, config["local_neighbors"]:]
nearest = np.linalg.norm(coords[local[:, 0]] - coords, axis=1)
scale = float(np.median(nearest[nearest > 0]))

adt = ad.read_h5ad(DATA_ROOT / "GBM/adt_dgat_clr.h5ad")[:, proteins]
target_truth = dense(adt.X)

# Compute target graph weights for calibration
x_mat = dense(rna.X)
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

# Load frozen calibration parameters for cross_holdout_GBM
with open(BASE / "SCP-NET-work/runs/protein_convex_rna_calibration_10_dataset_20261103/all_directions/cross_holdout_GBM/calibration_results.json", "r") as f:
    calib_cfg = json.load(f)

rho = calib_cfg["selected_rho"]
lambdas = np.asarray([calib_cfg["selected_lambda_by_protein"][p] for p in proteins])
target_weights = rho * sp_w + (1.0 - rho) * ft_w

# Load descriptors and PCA from run
ctx_npz = np.load(RUN_DIR / "rna_slice_context_descriptors.npz")
descriptors = ctx_npz["descriptors"]
canon_npz = np.load(RUN_DIR / "rna_slice_canonicalization.npz")
means = canon_npz["mean"]
stds = canon_npz["std"]

desc_t = torch.from_numpy(descriptors).to(device)
means_t = torch.from_numpy(means).to(device)
stds_t = torch.from_numpy(stds).to(device)

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

model = module.SliceConditionedNicheProtein(
    len(genes), len(proteins), config["hidden_dim"], config["heads"], config["ffn_dim"],
    config["dropout"], config["protein_embedding_dim"], config["router_rank"],
    context_input_dim=descriptors.shape[1], context_rank=config["context_rank"],
    protein_specific_router=not config["disable_protein_specific_router"],
    protein_aware_attention=not config["disable_protein_aware_attention"],
    use_slice_context=not config["disable_slice_context"],
    use_canonicalization=not config["disable_canonicalization"],
    use_context_adapter=not config["disable_context_adapter"],
    cognate_indices=cognate_indices, has_cognate=has_cognate,
    use_cognate_anchor=not config["disable_cognate_anchor"],
).to(device)

lp = relative_position(coords, local, scale)
rp = relative_position(coords, regional, scale)
target_idx = np.arange(len(x_mat), dtype=np.int64)
target_section_ids = np.full(len(target_idx), len(sources), dtype=np.int64)

results = []

# Load source validation history
hist_df = pd.read_csv(RUN_DIR / "metrics_history.csv")

for ep in range(1, 21):
    pth = CHECKPOINT_DIR / f"epoch_{ep:02d}.pth"
    state = torch.load(pth, map_location=device, weights_only=True)
    model.load_state_dict(state)
    model.eval()

    # Predict uncalibrated
    with torch.no_grad():
        preds = []
        for start in range(0, len(target_idx), config["batch_size"]):
            idx = target_idx[start:start + config["batch_size"]]
            inputs = module.model_batch_inputs(x_mat, local, regional, lp, rp, idx, target_section_ids, desc_t, means_t, stds_t, device)
            preds.append(model(*inputs, need_attention=False)["final"].cpu().numpy())
    pred_uncalib = np.concatenate(preds)

    # Apply calibration
    neighbor_pred = (target_weights[:, :, None] * pred_uncalib[n_nbrs]).sum(1)
    pred_calib = (1.0 - lambdas[None, :]) * pred_uncalib + lambdas[None, :] * neighbor_pred

    # Compute metrics
    m_uncalib = dgat_metrics(pred_uncalib, target_truth, proteins)
    m_calib = dgat_metrics(pred_calib, target_truth, proteins)

    src_val_pcc = hist_df.iloc[ep - 1]["source_validation_pcc_pro_mean"]
    src_val_rmse = hist_df.iloc[ep - 1]["source_validation_rmse_pro"]

    results.append({
        "epoch": ep,
        "src_val_PCC": src_val_pcc,
        "src_val_RMSE": src_val_rmse,
        "test_raw_PCC": m_uncalib["pcc_pro_mean"],
        "test_raw_RMSE": m_uncalib["rmse_global"],
        "test_calib_PCC": m_calib["pcc_pro_mean"],
        "test_calib_RMSE": m_calib["rmse_global"],
    })

res_df = pd.DataFrame(results)
res_df.to_csv(BASE / "SCP-NET-work/experiments/epoch_trajectory_probe/gbm_epoch_trajectory.csv", index=False)
print("=== Trajectory Evaluation Complete! ===")
