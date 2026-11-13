#!/usr/bin/env python3
"""Source-validation RNA-neighbour calibration on frozen ten-dataset models."""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import torch
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

BASE = Path(r"E:/DGAT")
DATA_ROOT = BASE / "ST2P/data/processed/DGAT_final10_common9384_muon017_axis0_20261101"
RUN_ROOT = BASE / "SCP-NET-work/runs/protein_convex_10_dataset_20261101"
PANEL_FILE = BASE / "SCP-NET-work/docs/ten_dataset_common_panel_20261101.json"
OUTPUT_ROOT = BASE / "SCP-NET-work/runs/protein_convex_rna_calibration_10_dataset_20261103"
EXPERIMENT_DIR = BASE / "ST2P/src"

sys.path.insert(0, str(EXPERIMENT_DIR))
from v9_robust_source_subset_test_oracle import corrected_knn, dense, panel
from triscale_center_protein_v2_relpos import relative_position
from triscale_center_protein_v9_dgat_compatible import dgat_metrics

RHO_GRID = [0.0, 0.25, 0.5, 0.75, 1.0]
LAMBDA_GRID = np.asarray([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], dtype=np.float64)
COGNATE = {
    "CD163": "CD163", "CR2": "CR2", "PCNA": "PCNA", "VIM": "VIM", "KRT5": "KRT5",
    "CD68": "CD68", "CEACAM8": "CEACAM8", "PTPRC_1": "PTPRC", "PTPRC_2": "PTPRC",
    "PAX5": "PAX5", "SDC1": "SDC1", "CD8A": "CD8A", "BCL2": "BCL2", "CD19": "CD19",
    "PDCD1": "PDCD1", "ACTA2": "ACTA2", "FCGR3A": "FCGR3A", "ITGAX": "ITGAX",
    "CXCR5": "CXCR5", "EPCAM": "EPCAM", "MS4A1": "MS4A1", "CD3E": "CD3E",
    "CD14": "CD14", "CD40": "CD40", "PECAM1": "PECAM1", "CD4": "CD4",
    "ITGAM": "ITGAM", "CD27": "CD27", "CCR7": "CCR7", "CD274": "CD274",
}


def load_snapshot(run_dir: Path):
    snapshot_path = run_dir / "script_snapshot.py"
    spec = importlib.util.spec_from_file_location(f"frozen_model_{run_dir.name}", snapshot_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load model snapshot: {snapshot_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def graph_weights(x: np.ndarray, coords: np.ndarray, seed: int):
    distances, neighbors = NearestNeighbors(n_neighbors=9).fit(coords).kneighbors(coords)
    distances, neighbors = distances[:, 1:], neighbors[:, 1:]
    n_components = min(32, x.shape[1], len(x) - 1)
    z = PCA(n_components=n_components, random_state=seed).fit_transform(x)
    feature_distances = np.linalg.norm(z[:, None, :] - z[neighbors], axis=-1)
    spatial_scale = np.median(distances[distances > 0])
    feature_scale = np.median(feature_distances[feature_distances > 0])
    spatial_weights = np.exp(-0.5 * (distances / (spatial_scale + 1e-8)) ** 2)
    feature_weights = np.exp(-0.5 * (feature_distances / (feature_scale + 1e-8)) ** 2)
    spatial_weights /= spatial_weights.sum(axis=1, keepdims=True)
    feature_weights /= feature_weights.sum(axis=1, keepdims=True)
    return neighbors, spatial_weights, feature_weights


def pcc_by_protein(pred: np.ndarray, truth: np.ndarray) -> np.ndarray:
    from scipy.stats import pearsonr
    vals = []
    for j in range(pred.shape[1]):
        a, b = pred[:, j], truth[:, j]
        vals.append(np.nan if np.std(a) <= 1e-8 or np.std(b) <= 1e-8 else float(pearsonr(a, b).statistic))
    return np.asarray(vals, dtype=np.float64)


def infer(model, base, x, local, regional, lp, rp, section_ids, descriptors, means, stds, indices, device, batch_size):
    model.eval()
    desc_t = torch.from_numpy(descriptors).to(device)
    means_t = torch.from_numpy(means).to(device)
    stds_t = torch.from_numpy(stds).to(device)
    output = []
    with torch.no_grad():
        for start in range(0, len(indices), batch_size):
            idx = np.asarray(indices[start:start + batch_size], dtype=np.int64)
            inputs = base.model_batch_inputs(x, local, regional, lp, rp, idx, section_ids, desc_t, means_t, stds_t, device)
            output.append(model(*inputs, need_attention=False)["final"].cpu().numpy())
    return np.concatenate(output).astype(np.float32, copy=False)


def prepare_and_calibrate(run_dir: Path, output_root: Path, exact_check: bool):
    result = json.loads((run_dir / "results.json").read_text(encoding="utf-8"))
    config = json.loads((run_dir / "config.json").read_text(encoding="utf-8"))
    sources, target_name = result["sources"], result["target"]
    seed = int(result["seed"])
    manifest = json.loads((DATA_ROOT / "prepared_manifest.json").read_text(encoding="utf-8"))
    dataset_names = [entry["dataset"] for entry in manifest["datasets"]]
    genes = panel(DATA_ROOT, dataset_names, "rna_dgat.h5ad")
    proteins = result["common_panel"]["proteins"]
    if len(genes) != 9384 or len(proteins) != 31:
        raise RuntimeError(f"wrong common panel: {len(genes)} genes, {len(proteins)} proteins")

    blocks = {}
    audit = {}
    for name in [*sources, target_name]:
        rna = ad.read_h5ad(DATA_ROOT / name / "rna_dgat.h5ad")[:, genes].copy()
        coords = np.asarray(rna.obsm["spatial"], dtype=np.float32)
        neighbors = corrected_knn(coords, int(config["regional_neighbors"]))
        local, regional = neighbors[:, :int(config["local_neighbors"])], neighbors[:, int(config["local_neighbors"]):]
        nearest = np.linalg.norm(coords[local[:, 0]] - coords, axis=1)
        scale = float(np.median(nearest[nearest > 0]))
        block = {
            "x": dense(rna.X), "local": local, "regional": regional,
            "lp": relative_position(coords, local, scale),
            "rp": relative_position(coords, regional, scale),
            "coords": coords, "obs": np.asarray(rna.obs_names),
        }
        if name in sources:
            adt = ad.read_h5ad(DATA_ROOT / name / "adt_dgat_clr.h5ad")[:, proteins].copy()
            if not rna.obs_names.equals(adt.obs_names):
                raise RuntimeError(f"RNA/ADT barcode mismatch: {name}")
            block["y"] = dense(adt.X)
        blocks[name] = block
        audit[name] = {"spots": len(block["x"]), "coordinate_scale": scale}

    split = np.load(run_dir / "source_split_indices.npz")
    fit_idx, val_idx = split["fit"], split["validation"]
    offsets = {}
    offset = 0
    parts = [[] for _ in range(6)]
    section_parts = []
    for pos, name in enumerate(sources):
        b = blocks[name]
        offsets[name] = (offset, offset + len(b["x"]))
        parts[0].append(b["x"]); parts[1].append(b["y"])
        parts[2].append(b["local"] + offset); parts[3].append(b["regional"] + offset)
        parts[4].append(b["lp"]); parts[5].append(b["rp"])
        section_parts.append(np.full(len(b["x"]), pos, dtype=np.int64))
        offset += len(b["x"])
    x, y, local, regional, lp, rp = [np.concatenate(p) for p in parts]
    section_ids = np.concatenate(section_parts)
    if int(fit_idx.max()) >= len(x) or int(val_idx.max()) >= len(x):
        raise RuntimeError(f"saved source split exceeds concatenated source spots: {run_dir.name}")

    # Recreate the RNA-only context features fitted by the original run on source-fit spots.
    context_components = min(int(config["context_pca_components"]), len(fit_idx) - 1, x.shape[1])
    context_pca = PCA(n_components=context_components, svd_solver="randomized", random_state=seed).fit(x[fit_idx])
    descriptors, means, stds = [], [], []
    for name in [*sources, target_name]:
        b = blocks[name]
        scores = context_pca.transform(b["x"])
        coord = b["lp"][:, 0, :2]
        summary = np.asarray([np.log1p(len(b["x"])), np.log1p(np.linalg.norm(coord, axis=1).mean()), np.log1p(np.linalg.norm(coord, axis=1).std())], dtype=np.float32)
        descriptors.append(np.concatenate((scores.mean(0), scores.std(0), summary)).astype(np.float32))
        means.append(b["x"].mean(0).astype(np.float32))
        stds.append(np.clip(b["x"].std(0), 0.1, None).astype(np.float32))
    descriptors = np.stack(descriptors)
    desc_center = descriptors[:-1].mean(0)
    desc_scale = np.clip(descriptors[:-1].std(0), 0.1, None)
    descriptors = np.clip((descriptors - desc_center) / desc_scale, -5.0, 5.0)
    means, stds = np.stack(means), np.stack(stds)

    cognate_map = {
        "CD163":"CD163", "CR2":"CR2", "PCNA":"PCNA", "VIM":"VIM", "KRT5":"KRT5", "CD68":"CD68", "CEACAM8":"CEACAM8",
        "PTPRC_1":"PTPRC", "PTPRC_2":"PTPRC", "PAX5":"PAX5", "SDC1":"SDC1", "CD8A":"CD8A", "BCL2":"BCL2", "CD19":"CD19",
        "PDCD1":"PDCD1", "ACTA2":"ACTA2", "FCGR3A":"FCGR3A", "ITGAX":"ITGAX", "CXCR5":"CXCR5", "EPCAM":"EPCAM", "MS4A1":"MS4A1",
        "CD3E":"CD3E", "CD14":"CD14", "CD40":"CD40", "PECAM1":"PECAM1", "CD4":"CD4", "ITGAM":"ITGAM", "CD27":"CD27", "CCR7":"CCR7", "CD274":"CD274",
    }
    gene_index = {g: i for i, g in enumerate(genes)}
    cognate_indices, has_cognate = [], []
    for protein in proteins:
        gene = cognate_map.get(protein)
        if gene in gene_index: cognate_indices.append(gene_index[gene]); has_cognate.append(1.0)
        else: cognate_indices.append(0); has_cognate.append(0.0)

    base = load_snapshot(run_dir)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = base.SliceConditionedNicheProtein(
        len(genes), len(proteins), int(config["hidden_dim"]), int(config["heads"]), int(config["ffn_dim"]),
        float(config["dropout"]), int(config["protein_embedding_dim"]), int(config["router_rank"]),
        context_input_dim=descriptors.shape[1], context_rank=int(config["context_rank"]),
        protein_specific_router=not config["disable_protein_specific_router"],
        protein_aware_attention=not config["disable_protein_aware_attention"],
        use_slice_context=not config["disable_slice_context"],
        use_canonicalization=not config["disable_canonicalization"],
        use_context_adapter=not config["disable_context_adapter"],
        cognate_indices=cognate_indices, has_cognate=has_cognate,
        use_cognate_anchor=not config["disable_cognate_anchor"],
    ).to(device)
    state = torch.load(run_dir / "best_source_validation.pth", map_location=device, weights_only=True)
    model.load_state_dict(state)

    batch_size = int(config["batch_size"])
    all_source_idx = np.arange(len(x), dtype=np.int64)
    source_pred = infer(model, base, x, local, regional, lp, rp, section_ids, descriptors, means, stds, all_source_idx, device, batch_size)
    target = blocks[target_name]
    target_idx = np.arange(len(target["x"]), dtype=np.int64)
    target_section = np.full(len(target_idx), len(sources), dtype=np.int64)
    target_pred = infer(model, base, target["x"], target["local"], target["regional"], target["lp"], target["rp"], target_section, descriptors, means, stds, target_idx, device, batch_size)

    saved_pred = np.load(run_dir / "target_source_selected_predictions.npy")
    max_reconstruction_diff = float(np.max(np.abs(target_pred - saved_pred)))
    if exact_check and max_reconstruction_diff > 2e-5:
        raise RuntimeError(f"base prediction reconstruction mismatch {run_dir.name}: max_abs={max_reconstruction_diff}")

    source_info = []
    for pos, name in enumerate(sources):
        lo, hi = offsets[name]
        local_val = val_idx[(val_idx >= lo) & (val_idx < hi)] - lo
        if len(local_val) == 0: raise RuntimeError(f"no source validation rows for {name}")
        idx, sw, fw = graph_weights(blocks[name]["x"], blocks[name]["coords"], seed)
        source_info.append((name, local_val.astype(np.int64), idx, sw, fw))

    score_details = []
    best_score, best_rho, best_lambdas = -np.inf, None, None
    for rho in RHO_GRID:
        scores = np.zeros((len(LAMBDA_GRID), len(proteins)), dtype=np.float64)
        counts = np.zeros_like(scores)
        for name, local_val, idx, sw, fw in source_info:
            lo, hi = offsets[name]
            pred = source_pred[lo:hi]
            truth = blocks[name]["y"]
            weights = rho * sw + (1.0 - rho) * fw
            neighbor_pred = (weights[:, :, None] * pred[idx]).sum(1)
            for li, lam in enumerate(LAMBDA_GRID):
                calibrated = (1.0 - lam) * pred + lam * neighbor_pred
                pcc = pcc_by_protein(calibrated[local_val], truth[local_val])
                valid = np.isfinite(pcc)
                scores[li, valid] += pcc[valid]
                counts[li, valid] += 1
        mean_scores = scores / np.maximum(counts, 1)
        chosen = np.argmax(mean_scores, axis=0)
        selected = mean_scores[chosen, np.arange(len(proteins))]
        value = float(np.nanmean(selected))
        score_details.append({"rho": rho, "source_validation_macro_pcc": value})
        if value > best_score:
            best_score = value
            best_rho = rho
            best_lambdas = LAMBDA_GRID[chosen]

    target_idx_graph, target_sw, target_fw = graph_weights(target["x"], target["coords"], seed)
    target_weights = best_rho * target_sw + (1.0 - best_rho) * target_fw
    target_neighbor_pred = (target_weights[:, :, None] * target_pred[target_idx_graph]).sum(1)
    calibrated_pred = ((1.0 - best_lambdas[None, :]) * target_pred + best_lambdas[None, :] * target_neighbor_pred).astype(np.float32)

    # Target protein values are opened only after all calibration choices are frozen.
    adt = ad.read_h5ad(DATA_ROOT / target_name / "adt_dgat_clr.h5ad")[:, proteins].copy()
    target_obs = np.asarray(target["obs"]).astype(str)
    if not np.array_equal(np.asarray(adt.obs_names).astype(str), target_obs):
        raise RuntimeError(f"target RNA/ADT barcode order mismatch: {target_name}")
    target_truth = dense(adt.X)
    base_metrics = dgat_metrics(target_pred, target_truth, proteins)
    calibrated_metrics = dgat_metrics(calibrated_pred, target_truth, proteins)

    out_dir = output_root / run_dir.name
    out_dir.mkdir(parents=True, exist_ok=False)
    np.save(out_dir / "rna_calibrated_predictions.npy", calibrated_pred)
    pd = __import__("pandas")
    pd.DataFrame(calibrated_metrics["per_protein"]).to_csv(out_dir / "rna_calibrated_per_protein_metrics.csv", index=False)
    record = {
        "direction": run_dir.name, "target": target_name, "sources": sources,
        "checkpoint": str(run_dir / "best_source_validation.pth"), "best_epoch": result["best_epoch"],
        "source_validation_only_selection": True, "target_adt_opened_after_selection": True,
        "common_panel": {"genes": len(genes), "proteins": proteins},
        "rho_grid": RHO_GRID, "lambda_grid": LAMBDA_GRID.tolist(),
        "selected_rho": float(best_rho), "selected_lambda_by_protein": {p: float(v) for p, v in zip(proteins, best_lambdas)},
        "source_validation_macro_pcc_by_rho": score_details,
        "baseline_reconstruction_max_abs_diff": max_reconstruction_diff,
        "baseline_metrics": {k: v for k, v in base_metrics.items() if k != "per_protein"},
        "calibrated_metrics": {k: v for k, v in calibrated_metrics.items() if k != "per_protein"},
        "baseline_rmse_pro": float(np.nanmean([r["rmse"] for r in base_metrics["per_protein"]])),
        "baseline_mae_pro": float(np.nanmean([r["mae"] for r in base_metrics["per_protein"]])),
        "calibrated_rmse_pro": float(np.nanmean([r["rmse"] for r in calibrated_metrics["per_protein"]])),
        "calibrated_mae_pro": float(np.nanmean([r["mae"] for r in calibrated_metrics["per_protein"]])),
        "prediction_file": "rna_calibrated_predictions.npy",
        "truth_file": str(run_dir / "target_source_selected_truth.npy"),
    }
    (out_dir / "calibration_results.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--direction", default="all", help="all or one exact ProteinConvex result directory name")
    parser.add_argument("--exact-check", action="store_true")
    parser.add_argument("--run-root", type=Path, default=RUN_ROOT, help="directory containing frozen direction runs")
    parser.add_argument("--output-root", type=Path, default=None, help="independent output directory; defaults to the baseline calibration root")
    args = parser.parse_args()
    if args.output_root is None:
        output_root = OUTPUT_ROOT / ("pilot" if args.direction != "all" else "all_directions")
    else:
        output_root = args.output_root
    if output_root.exists(): raise FileExistsError(f"refusing to overwrite {output_root}")
    output_root.mkdir(parents=True)
    run_dirs = sorted(p for p in args.run_root.iterdir() if p.is_dir() and (p / "results.json").exists())
    if args.direction != "all": run_dirs = [p for p in run_dirs if p.name == args.direction]
    if not run_dirs: raise RuntimeError(f"no run matched {args.direction}")
    records = []
    for run_dir in run_dirs:
        print(f"START {run_dir.name}", flush=True)
        record = prepare_and_calibrate(run_dir, output_root, args.exact_check)
        records.append(record)
        print(json.dumps({"direction": record["direction"], "target": record["target"], "rho": record["selected_rho"], "lambda_mean": float(np.mean(list(record["selected_lambda_by_protein"].values()))), "source_val": max(r["source_validation_macro_pcc"] for r in record["source_validation_macro_pcc_by_rho"]), "baseline": record["baseline_metrics"], "calibrated": record["calibrated_metrics"], "prediction_diff": record["baseline_reconstruction_max_abs_diff"]}, ensure_ascii=False), flush=True)
    summary = {
        "candidate": "ProteinConvex + source-only RNA neighbor calibration",
        "panel": {"genes": 9384, "proteins": 31},
        "selection_protocol": "per-direction source validation only; target ADT loaded after rho/lambda selection",
        "direction_count": len(records), "records": records,
    }
    (output_root / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"COMPLETE directions={len(records)} output={output_root}", flush=True)


if __name__ == "__main__":
    main()
