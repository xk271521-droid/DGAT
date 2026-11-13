#!/usr/bin/env python3
"""V9 under the official DGAT RNA/ADT preprocessing protocol.

This is a deliberately separate benchmark entrypoint.  It consumes raw paired
RNA/ADT H5ADs, applies the DGAT CytAssist QC and transforms once, then trains
the protected V9 architecture.  It never modifies legacy V9 artefacts.
"""
from __future__ import annotations

import argparse
import json
import random
from copy import deepcopy
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from scipy.stats import pearsonr, spearmanr
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

from triscale_center_protein_v2_relpos import batch_inputs, compute_metrics, dense, relative_position, row_minmax
from triscale_center_protein_v5_program_decoder import ProgramResidualDecoder
from triscale_center_protein_v7_anchored_cross_attention import NeighbourCrossPool


def clean_protein_names(var_names: pd.Index) -> list[str]:
    """Match DGAT's public name cleaning before panel intersection/QC."""
    cleaned = [str(v).replace("-", "_") for v in var_names]
    groups: dict[str, list[str]] = {}
    parsed: list[tuple[str, str | None, str]] = []
    for name in cleaned:
        base, suffix = name.rsplit("_", 1) if "_" in name and name.rsplit("_", 1)[1].isdigit() else (name, None)
        groups.setdefault(base, []).append(name)
        parsed.append((base, suffix, name))
    return [base if suffix is not None and len(groups[base]) == 1 else original
            for base, suffix, original in parsed]


def harmonize_adt_feature_identifiers(adt: ad.AnnData) -> ad.AnnData:
    """Resolve vendor aliases before intersecting ADT panels across sections.

    The CytAssist matrices used for LN, Tonsil and Breast label the two CD45
    antibodies with the duplicate display name ``PTPRC`` while retaining their
    unique vendor feature identifiers in either ``var['gene_ids']`` or
    ``var['feature_id']`` (``PTPRC_1`` and ``PTPRC_2``).  Human Tonsil AddOns
    exposes those identifiers as its display names.  Intersecting display names
    alone would consequently discard both shared antibody channels.  Use the
    vendor identifier only for these two documented aliases, retain every other
    display name, and preserve the original label in ``var`` for auditability.
    """
    adt = adt.copy()
    original = adt.var_names.astype(str).copy()
    adt.var["panel_display_name_original"] = original.to_numpy()
    identifier_column = next(
        (column for column in ("gene_ids", "feature_id") if column in adt.var.columns),
        None,
    )
    if identifier_column is None:
        return adt
    vendor_ids = adt.var[identifier_column].astype(str).to_numpy()
    harmonized = [
        vendor if vendor in {"PTPRC_1", "PTPRC_2"} else display
        for display, vendor in zip(original, vendor_ids)
    ]
    adt.var["panel_feature_identifier"] = harmonized
    adt.var["panel_feature_identifier_source"] = identifier_column
    adt.var_names = pd.Index(harmonized)
    return adt


def dgat_preprocess(rna: ad.AnnData, adt: ad.AnnData) -> tuple[ad.AnnData, ad.AnnData]:
    """Apply DGAT's public CytAssist QC, RNA transform and ADT CLR once.

    The global RNA scaling statistics use all spots, as in DGAT's own
    ``Preprocessing.normalize``.  This is transductive only with respect to
    RNA input values; protein labels are not used for the transform.
    """
    shared = rna.obs_names.intersection(adt.obs_names, sort=False)
    rna, adt = rna[shared].copy(), adt[shared].copy()
    if "spatial" not in rna.obsm:
        raise KeyError("RNA H5AD must contain obsm['spatial']")
    if np.min(dense(rna.X)) < 0 or np.min(dense(adt.X)) < 0:
        raise ValueError("DGAT preprocessing expects raw non-negative count matrices.")

    # The raw 10x-derived registry preserves vendor gene symbols, which may
    # repeat.  AnnData requires a unique axis for the DGAT QC subset; this only
    # disambiguates duplicate names and neither drops nor changes count values.
    rna.var_names_make_unique()
    adt = harmonize_adt_feature_identifiers(adt)
    adt.var_names = clean_protein_names(adt.var_names)
    adt.var_names_make_unique()
    rna.var["mt"] = rna.var_names.astype(str).str.startswith("MT-")
    sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], inplace=True, percent_top=None)
    original_rna = rna.copy()
    min_cells = max(1, int(np.ceil(rna.n_obs * 0.025)))
    sc.pp.filter_genes(rna, min_cells=min_cells)
    # DGAT protects genes cognate to the retained protein panel from its
    # expression-frequency filter.
    gene_set = set(rna.var_names.astype(str))
    common_prefixes = {str(p).split("_")[0] for p in adt.var_names}
    protected = {g for g in original_rna.var_names.astype(str) if g.split("_")[0] in common_prefixes}
    retained = [g for g in original_rna.var_names.astype(str) if g in gene_set or g in protected]
    rna = original_rna[:, retained].copy()
    sc.pp.filter_cells(rna, min_genes=700)
    rna = rna[rna.obs["pct_counts_mt"] < 35.0].copy()
    adt = adt[rna.obs_names].copy()
    isotype = (adt.var_names.astype(str).str.startswith("mouse_") |
               adt.var_names.astype(str).str.startswith("rat_") |
               adt.var_names.astype(str).str.startswith("mouse.") |
               adt.var_names.astype(str).str.startswith("rat."))
    adt = adt[:, ~isotype].copy()
    rna.layers["dgat_raw"] = rna.X.copy()
    adt.layers["dgat_raw"] = adt.X.copy()
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    sc.pp.scale(rna, max_value=10)
    # ``muon.prot.pp.clr`` used by DGAT is log1p counts centered within spot.
    protein = np.log1p(dense(adt.X).astype(np.float32))
    adt.X = (protein - protein.mean(axis=1, keepdims=True)).astype(np.float32)
    return rna, adt


def dgat_metrics(prediction: np.ndarray, target: np.ndarray, proteins: list[str]) -> dict:
    """DGAT's correlation/RMSE definitions plus explicit per-protein values."""
    protein_rows, spot_pcc, spot_spearman = [], [], []
    for col, protein in enumerate(proteins):
        p, t = prediction[:, col], target[:, col]
        pearson = np.nan if np.std(p) == 0 or np.std(t) == 0 else float(pearsonr(p, t)[0])
        spearman = np.nan if np.std(p) == 0 or np.std(t) == 0 else float(spearmanr(p, t)[0])
        protein_rows.append({"protein": protein, "pcc": pearson, "spearman": spearman,
                             "rmse": float(np.sqrt(np.mean((p - t) ** 2))),
                             "mae": float(np.mean(np.abs(p - t)))})
    for p, t in zip(prediction, target):
        if np.std(p) and np.std(t):
            spot_pcc.append(float(pearsonr(p, t)[0]))
            spot_spearman.append(float(spearmanr(p, t)[0]))
    return {
        "pcc_pro_mean": float(np.nanmean([r["pcc"] for r in protein_rows])),
        "pcc_spot_mean": float(np.nanmean(spot_pcc)),
        "spearman_pro_mean": float(np.nanmean([r["spearman"] for r in protein_rows])),
        "spearman_spot_mean": float(np.nanmean(spot_spearman)),
        "rmse_global": float(np.sqrt(np.mean((prediction - target) ** 2))),
        "mae_global": float(np.mean(np.abs(prediction - target))),
        "per_protein": protein_rows,
    }


def evaluate_dgat(model, x, y, local_index, regional_index, local_position, regional_position,
                  index, proteins, device, batch_size, save_details=False):
    model.eval(); predictions, activities, residuals, local_attn, regional_attn, scale_attn = [], [], [], [], [], []
    with torch.no_grad():
        for start in range(0, len(index), batch_size):
            idx = index[start:start + batch_size]
            out = model(*batch_inputs(x, local_index, regional_index, local_position, regional_position, idx, device),
                        need_attention=save_details)
            predictions.append(out["final"].cpu().numpy())
            if save_details:
                activities.append(out["program_activity"].cpu().numpy()); residuals.append(out["residual"].cpu().numpy())
                local_attn.append(out["local_attention"].cpu().numpy()); regional_attn.append(out["regional_attention"].cpu().numpy())
                scale_attn.append(out["scale_attention"].cpu().numpy())
    prediction = np.concatenate(predictions)
    details = None if not save_details else {"program_activity": np.concatenate(activities), "residual": np.concatenate(residuals),
        "local_attention": np.concatenate(local_attn), "regional_attention": np.concatenate(regional_attn),
        "scale_attention": np.concatenate(scale_attn)}
    return dgat_metrics(prediction, y[index], proteins), prediction, details


class TriScaleCenterProteinV9(nn.Module):
    """No center token can enter neighbour values; scale routing remains explicit."""

    def __init__(self, n_genes: int, basis: np.ndarray, mean: np.ndarray, hidden: int,
                 heads: int, ffn_dim: int, dropout: float):
        super().__init__()
        if hidden % heads:
            raise ValueError("hidden_dim must be divisible by heads")
        self.rna_encoder = nn.Sequential(nn.Linear(n_genes, hidden), nn.ReLU(), nn.Dropout(dropout))
        self.center_norm = nn.LayerNorm(hidden)
        self.position_encoder = nn.Sequential(
            nn.Linear(3, hidden), nn.ReLU(), nn.Linear(hidden, hidden, bias=False),
        )
        self.local_pool = NeighbourCrossPool(hidden, heads, ffn_dim, dropout)
        self.regional_pool = NeighbourCrossPool(hidden, heads, ffn_dim, dropout)
        # Same minimal router used successfully in V5: [Center, Local, Regional] compete directly.
        self.level_router = nn.Sequential(nn.Linear(hidden * 3, 3), nn.Softmax(dim=1))
        self.protein_program_decoder = ProgramResidualDecoder(hidden, basis, mean, dropout)

    def encode_neighbours(self, x: torch.Tensor, position: torch.Tensor) -> torch.Tensor:
        batch, token_count, genes = x.shape
        rna = self.rna_encoder(x.reshape(batch * token_count, genes)).reshape(batch, token_count, -1)
        return rna + self.position_encoder(position)

    def forward(self, center, local, regional, local_position, regional_position, need_attention: bool = False):
        z_center = self.center_norm(self.rna_encoder(center))
        query = z_center.unsqueeze(1)
        local_tokens = self.encode_neighbours(local, local_position)
        regional_tokens = self.encode_neighbours(regional, regional_position)
        z_local, local_attention = self.local_pool(query, local_tokens, need_attention)
        z_regional, regional_attention = self.regional_pool(query, regional_tokens, need_attention)
        z_local, z_regional = z_local.squeeze(1), z_regional.squeeze(1)
        levels = torch.stack((z_center, z_local, z_regional), dim=1)
        level_weights = self.level_router(torch.cat((z_center, z_local, z_regional), dim=1))
        z_final = (levels * level_weights.unsqueeze(-1)).sum(dim=1)
        final, activity, residual = self.protein_program_decoder(z_final)
        center_pred, _, _ = self.protein_program_decoder(z_center)
        local_pred, _, _ = self.protein_program_decoder(z_local)
        regional_pred, _, _ = self.protein_program_decoder(z_regional)
        # The V7 evaluator stores this common key; add singleton axes for a uniform artifact schema.
        return {
            "final": final, "center": center_pred, "local": local_pred, "regional": regional_pred,
            "program_activity": activity, "residual": residual,
            "local_attention": local_attention, "regional_attention": regional_attention,
            "scale_attention": level_weights[:, None, None, :],
        }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna_h5ad", required=True); ap.add_argument("--adt_h5ad", required=True); ap.add_argument("--out_dir", required=True)
    ap.add_argument("--dataset_name", default="Human Lymph Node A1"); ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--train_ratio", type=float, default=0.9); ap.add_argument("--validation_ratio", type=float, default=0.1)
    ap.add_argument("--local_neighbors", type=int, default=6); ap.add_argument("--regional_neighbors", type=int, default=24)
    ap.add_argument("--epochs", type=int, default=500); ap.add_argument("--patience", type=int, default=50)
    ap.add_argument("--batch_size", type=int, default=64); ap.add_argument("--hidden_dim", type=int, default=256)
    ap.add_argument("--heads", type=int, default=4); ap.add_argument("--ffn_dim", type=int, default=256); ap.add_argument("--dropout", type=float, default=0.3)
    ap.add_argument("--lr", type=float, default=5e-4); ap.add_argument("--weight_decay", type=float, default=1e-4)
    ap.add_argument("--aux_weight", type=float, default=0.25); ap.add_argument("--program_components", type=int, default=9)
    args = ap.parse_args()
    random.seed(args.seed); np.random.seed(args.seed); torch.manual_seed(args.seed)
    if torch.cuda.is_available(): torch.cuda.manual_seed_all(args.seed)
    torch.backends.cudnn.deterministic = True; torch.backends.cudnn.benchmark = False

    out_dir = Path(args.out_dir)
    if out_dir.exists(): raise FileExistsError(f"Refusing to overwrite existing result directory: {out_dir}")
    out_dir.mkdir(parents=True); (out_dir / "config.json").write_text(json.dumps(vars(args), indent=2), encoding="utf-8")
    raw_rna, raw_adt = ad.read_h5ad(args.rna_h5ad), ad.read_h5ad(args.adt_h5ad)
    rna, adt = dgat_preprocess(raw_rna, raw_adt)
    # The RNA and ADT transforms have already occurred above.  In particular,
    # do not apply V9's legacy row-min--max transform after DGAT scaling.
    x, y = dense(rna.X).astype(np.float32), dense(adt.X).astype(np.float32); keep = np.var(y, axis=0) > 0
    y, proteins = y[:, keep], [str(p) for p in adt.var_names[keep]]; coords = np.asarray(rna.obsm["spatial"], dtype=np.float32)
    (out_dir / "preprocessing_audit.json").write_text(json.dumps({
        "rna_input": "DGAT: QC (>=700 genes, mt<35%, genes in >=2.5% spots plus cognate genes), normalize_total(1e4), log1p, scale, clip=10",
        "protein_target": "DGAT CLR: log1p(raw ADT) - mean(log1p(raw ADT)) within each spot; mouse_/rat_ controls removed",
        "raw_rna_shape": [int(raw_rna.n_obs), int(raw_rna.n_vars)],
        "raw_adt_shape": [int(raw_adt.n_obs), int(raw_adt.n_vars)],
        "post_qc_rna_shape": [int(rna.n_obs), int(rna.n_vars)],
        "post_qc_adt_shape": [int(adt.n_obs), int(adt.n_vars)],
        "target_shape": [int(y.shape[0]), int(y.shape[1])],
        "target_min": float(y.min()), "target_max": float(y.max()),
    }, indent=2), encoding="utf-8")
    all_neighbors = NearestNeighbors(n_neighbors=args.regional_neighbors + 1).fit(coords).kneighbors(return_distance=False)[:, 1:]
    local_index, regional_index = all_neighbors[:, :args.local_neighbors], all_neighbors[:, args.local_neighbors:args.regional_neighbors]
    nearest_distance = np.linalg.norm(coords[local_index[:, 0]] - coords, axis=1); coordinate_scale = float(np.median(nearest_distance[nearest_distance > 0]))
    local_position = relative_position(coords, local_index, coordinate_scale); regional_position = relative_position(coords, regional_index, coordinate_scale)

    rng = np.random.default_rng(args.seed); shuffled = np.arange(len(x)); rng.shuffle(shuffled); outer_cut = int(len(shuffled) * args.train_ratio)
    outer_train, test_index = np.sort(shuffled[:outer_cut]), np.sort(shuffled[outer_cut:]); inner = outer_train.copy(); rng = np.random.default_rng(args.seed); rng.shuffle(inner)
    val_n = max(1, int(round(len(inner) * args.validation_ratio))); validation_index, fit_index = np.sort(inner[:val_n]), np.sort(inner[val_n:])
    np.savez_compressed(out_dir / "split_indices.npz", fit=fit_index, validation=validation_index, test=test_index)
    n_components = min(args.program_components, len(fit_index), y.shape[1]); pca = PCA(n_components=n_components, svd_solver="full", random_state=args.seed).fit(y[fit_index])
    np.save(out_dir / "pca_basis_fit.npy", pca.components_.astype(np.float32)); np.save(out_dir / "pca_mean_fit.npy", pca.mean_.astype(np.float32))
    (out_dir / "pca_fit_audit.json").write_text(json.dumps({"fit_spots": int(len(fit_index)), "components": int(n_components), "explained_variance_ratio": pca.explained_variance_ratio_.tolist(), "cumulative_explained_variance": float(pca.explained_variance_ratio_.sum()), "basis_fitted_on": "fit proteins only"}, indent=2), encoding="utf-8")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = TriScaleCenterProteinV9(x.shape[1], pca.components_.astype(np.float32), pca.mean_.astype(np.float32), args.hidden_dim, args.heads, args.ffn_dim, args.dropout).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay); loader = DataLoader(TensorDataset(torch.from_numpy(fit_index)), batch_size=args.batch_size, shuffle=True, generator=torch.Generator().manual_seed(args.seed))
    print(f"TriScaleCenterProtein_V9_DGATCompatible parameters: {sum(p.numel() for p in model.parameters()):,}; pca_variance={pca.explained_variance_ratio_.sum():.4f}; coordinate_scale={coordinate_scale:.4f}", flush=True)
    mse = nn.MSELoss(); best_pcc, best_epoch, stale, best_state, history = -np.inf, 0, 0, None, []
    for epoch in range(1, args.epochs + 1):
        model.train(); losses = []
        for (batch_index,) in loader:
            idx = batch_index.numpy(); target = torch.from_numpy(y[idx]).to(device, non_blocking=True)
            out = model(*batch_inputs(x, local_index, regional_index, local_position, regional_position, idx, device))
            auxiliary = (mse(out["center"], target) + mse(out["local"], target) + mse(out["regional"], target)) / 3.0
            loss = mse(out["final"], target) + args.aux_weight * auxiliary
            optimizer.zero_grad(set_to_none=True); loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step(); losses.append(float(loss.detach().cpu()))
        val_metric, _, _ = evaluate_dgat(model, x, y, local_index, regional_index, local_position, regional_position, validation_index, proteins, device, args.batch_size)
        improved = val_metric["pcc_pro_mean"] > best_pcc + 1e-6
        if improved:
            best_pcc, best_epoch, stale = val_metric["pcc_pro_mean"], epoch, 0; best_state = deepcopy({k: v.detach().cpu() for k, v in model.state_dict().items()}); torch.save(best_state, out_dir / "best_validation.pth")
        else: stale += 1
        history.append({"epoch": epoch, "train_loss": float(np.mean(losses)), "validation_pcc_pro_mean": val_metric["pcc_pro_mean"], "validation_rmse_global": val_metric["rmse_global"], "validation_mae_global": val_metric["mae_global"], "is_best": improved})
        pd.DataFrame(history).to_csv(out_dir / "metrics_history.csv", index=False)
        print(f"Epoch [{epoch:03d}/{args.epochs}]" + (" [BEST]" if improved else "") + f" loss={np.mean(losses):.5f} validation_PCCPro={val_metric['pcc_pro_mean']:.4f}", flush=True)
        if stale >= args.patience: print(f"Early stopped after {args.patience} epochs without validation PCC improvement.", flush=True); break
    model.load_state_dict(best_state)
    test_metric, prediction, details = evaluate_dgat(model, x, y, local_index, regional_index, local_position, regional_position, test_index, proteins, device, args.batch_size, save_details=True)
    np.save(out_dir / "test_predictions.npy", prediction); np.savez_compressed(out_dir / "test_v9_details.npz", test_index=test_index, **details); pd.DataFrame(test_metric.pop("per_protein")).to_csv(out_dir / "per_protein_test_metrics.csv", index=False)
    result = {"model": "TriScaleCenterProtein_V9_DGATCompatible", "setting": "same-slice strict spot holdout; whole-slice transductive RNA context and DGAT global RNA scaling", "selection_protocol": "train-side validation mean per-protein Pearson only", "rna_preprocessing": "DGAT QC + total-count normalization + log1p + gene-wise scale + clip10", "protein_target_preprocessing": "DGAT CLR", "pca_protocol": "fixed basis and mean fitted on fit proteins only", "dataset": args.dataset_name, "seed": args.seed, "n_spots": int(len(x)), "n_genes": int(x.shape[1]), "n_proteins": int(y.shape[1]), "n_fit": int(len(fit_index)), "n_validation": int(len(validation_index)), "n_test": int(len(test_index)), "local_neighbors": int(local_index.shape[1]), "regional_outer_neighbors": int(regional_index.shape[1]), "coordinate_scale": coordinate_scale, "program_components": int(n_components), "pca_explained_variance": float(pca.explained_variance_ratio_.sum()), "n_parameters": int(sum(p.numel() for p in model.parameters())), "best_epoch": int(best_epoch), "best_validation_pcc_pro_mean": float(best_pcc), "test": test_metric}
    (out_dir / "results.json").write_text(json.dumps(result, indent=2), encoding="utf-8"); print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
