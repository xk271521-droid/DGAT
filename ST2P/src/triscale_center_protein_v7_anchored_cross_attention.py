#!/usr/bin/env python3
"""V7: neighbour-only cross-attention with a center-anchored scale decoder.

V6 showed that placing the center token in every self-attention key/value set
permits a center-copy shortcut.  V7 keeps the center only as a query: local
and regional summaries are necessarily formed from real neighbour RNA tokens.
The final fusion retains the center identity through a residual path, while its
cross-attention keys/values are the two neighbour-derived scale summaries.
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
import torch
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

from triscale_center_protein_v2_relpos import batch_inputs, compute_metrics, dense, relative_position, row_minmax
from triscale_center_protein_v5_program_decoder import ProgramResidualDecoder


class NeighbourCrossPool(nn.Module):
    """Center-query attention over neighbour-only keys and values.

    The returned representation has no center residual.  This prevents the
    local/regional auxiliary branches from satisfying their objective by simply
    copying the center RNA representation.
    """

    def __init__(self, hidden: int, heads: int, ffn_dim: int, dropout: float):
        super().__init__()
        self.query_norm = nn.LayerNorm(hidden)
        self.context_norm = nn.LayerNorm(hidden)
        self.attention = nn.MultiheadAttention(hidden, heads, dropout=dropout, batch_first=True)
        self.norm_attention = nn.LayerNorm(hidden)
        self.ffn = nn.Sequential(
            nn.Linear(hidden, ffn_dim), nn.GELU(), nn.Dropout(dropout), nn.Linear(ffn_dim, hidden),
        )
        self.norm_ffn = nn.LayerNorm(hidden)
        self.dropout = nn.Dropout(dropout)

    def forward(self, query: torch.Tensor, context: torch.Tensor, need_attention: bool = False):
        message, weight = self.attention(
            self.query_norm(query), self.context_norm(context), self.context_norm(context),
            need_weights=need_attention, average_attn_weights=False,
        )
        output = self.norm_attention(message)
        output = self.norm_ffn(output + self.dropout(self.ffn(output)))
        return output, weight


class ScaleCrossFusion(nn.Module):
    """Fuse neighbour-derived local/regional evidence into the center anchor."""

    def __init__(self, hidden: int, heads: int, ffn_dim: int, dropout: float):
        super().__init__()
        self.query_norm = nn.LayerNorm(hidden)
        self.context_norm = nn.LayerNorm(hidden)
        self.attention = nn.MultiheadAttention(hidden, heads, dropout=dropout, batch_first=True)
        self.norm_attention = nn.LayerNorm(hidden)
        self.ffn = nn.Sequential(
            nn.Linear(hidden, ffn_dim), nn.GELU(), nn.Dropout(dropout), nn.Linear(ffn_dim, hidden),
        )
        self.norm_ffn = nn.LayerNorm(hidden)
        self.dropout = nn.Dropout(dropout)

    def forward(self, center: torch.Tensor, scale_context: torch.Tensor, need_attention: bool = False):
        message, weight = self.attention(
            self.query_norm(center), self.context_norm(scale_context), self.context_norm(scale_context),
            need_weights=need_attention, average_attn_weights=False,
        )
        output = self.norm_attention(center + self.dropout(message))
        output = self.norm_ffn(output + self.dropout(self.ffn(output)))
        return output, weight


class TriScaleCenterProteinV7(nn.Module):
    """V5 decoder with neighbour-only intra-scale and anchored inter-scale attention."""

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
        self.scale_fusion = ScaleCrossFusion(hidden, heads, ffn_dim, dropout)
        self.protein_program_decoder = ProgramResidualDecoder(hidden, basis, mean, dropout)

    def encode_neighbours(self, x: torch.Tensor, position: torch.Tensor) -> torch.Tensor:
        batch, token_count, genes = x.shape
        rna = self.rna_encoder(x.reshape(batch * token_count, genes)).reshape(batch, token_count, -1)
        return rna + self.position_encoder(position)

    def forward(self, center, local, regional, local_position, regional_position, need_attention: bool = False):
        z_center = self.center_norm(self.rna_encoder(center))                 # [B, 256]
        query = z_center.unsqueeze(1)                                         # [B, 1, 256]
        local_tokens = self.encode_neighbours(local, local_position)          # [B, 6, 256]
        regional_tokens = self.encode_neighbours(regional, regional_position) # [B, 18, 256]

        z_local, local_attention = self.local_pool(query, local_tokens, need_attention)
        z_regional, regional_attention = self.regional_pool(query, regional_tokens, need_attention)
        z_local, z_regional = z_local.squeeze(1), z_regional.squeeze(1)       # neighbour-derived only
        z_final, scale_attention = self.scale_fusion(
            query, torch.stack((z_local, z_regional), dim=1), need_attention,
        )
        z_final = z_final.squeeze(1)                                          # [B, 256]

        final, activity, residual = self.protein_program_decoder(z_final)
        center_pred, _, _ = self.protein_program_decoder(z_center)
        local_pred, _, _ = self.protein_program_decoder(z_local)
        regional_pred, _, _ = self.protein_program_decoder(z_regional)
        return {
            "final": final, "center": center_pred, "local": local_pred, "regional": regional_pred,
            "program_activity": activity, "residual": residual,
            "local_attention": local_attention, "regional_attention": regional_attention,
            "scale_attention": scale_attention,
        }


def evaluate(model, x, y, local_index, regional_index, local_position, regional_position,
             index, proteins, device, batch_size, save_details=False):
    model.eval()
    predictions, activities, residuals, local_attn, regional_attn, scale_attn = [], [], [], [], [], []
    with torch.no_grad():
        for start in range(0, len(index), batch_size):
            idx = index[start:start + batch_size]
            out = model(*batch_inputs(x, local_index, regional_index, local_position, regional_position, idx, device),
                        need_attention=save_details)
            predictions.append(out["final"].cpu().numpy())
            if save_details:
                activities.append(out["program_activity"].cpu().numpy())
                residuals.append(out["residual"].cpu().numpy())
                local_attn.append(out["local_attention"].cpu().numpy())
                regional_attn.append(out["regional_attention"].cpu().numpy())
                scale_attn.append(out["scale_attention"].cpu().numpy())
    prediction = np.concatenate(predictions)
    details = None
    if save_details:
        details = {
            "program_activity": np.concatenate(activities), "residual": np.concatenate(residuals),
            "local_attention": np.concatenate(local_attn), "regional_attention": np.concatenate(regional_attn),
            "scale_attention": np.concatenate(scale_attn),
        }
    return compute_metrics(prediction, y[index], proteins), prediction, details


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
    rna, adt = ad.read_h5ad(args.rna_h5ad), ad.read_h5ad(args.adt_h5ad)
    shared = rna.obs_names.intersection(adt.obs_names, sort=False); rna, adt = rna[shared].copy(), adt[shared].copy()
    if "spatial" not in rna.obsm: raise KeyError("RNA h5ad must contain obsm['spatial']")
    x, y = row_minmax(dense(rna.X)), row_minmax(dense(adt.X)); keep = np.var(y, axis=0) > 0
    y, proteins = y[:, keep], [str(p) for p in adt.var_names[keep]]; coords = np.asarray(rna.obsm["spatial"], dtype=np.float32)
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
    model = TriScaleCenterProteinV7(x.shape[1], pca.components_.astype(np.float32), pca.mean_.astype(np.float32), args.hidden_dim, args.heads, args.ffn_dim, args.dropout).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay); loader = DataLoader(TensorDataset(torch.from_numpy(fit_index)), batch_size=args.batch_size, shuffle=True, generator=torch.Generator().manual_seed(args.seed))
    print(f"TriScaleCenterProtein_V7_AnchoredCrossAttention parameters: {sum(p.numel() for p in model.parameters()):,}; pca_variance={pca.explained_variance_ratio_.sum():.4f}; coordinate_scale={coordinate_scale:.4f}", flush=True)
    mse = nn.MSELoss(); best_pcc, best_epoch, stale, best_state, history = -np.inf, 0, 0, None, []
    for epoch in range(1, args.epochs + 1):
        model.train(); losses = []
        for (batch_index,) in loader:
            idx = batch_index.numpy(); target = torch.from_numpy(y[idx]).to(device, non_blocking=True)
            out = model(*batch_inputs(x, local_index, regional_index, local_position, regional_position, idx, device))
            auxiliary = (mse(out["center"], target) + mse(out["local"], target) + mse(out["regional"], target)) / 3.0
            loss = mse(out["final"], target) + args.aux_weight * auxiliary
            optimizer.zero_grad(set_to_none=True); loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step(); losses.append(float(loss.detach().cpu()))
        val_metric, _, _ = evaluate(model, x, y, local_index, regional_index, local_position, regional_position, validation_index, proteins, device, args.batch_size)
        improved = val_metric["pcc"] > best_pcc + 1e-6
        if improved:
            best_pcc, best_epoch, stale = val_metric["pcc"], epoch, 0; best_state = deepcopy({k: v.detach().cpu() for k, v in model.state_dict().items()}); torch.save(best_state, out_dir / "best_validation.pth")
        else: stale += 1
        history.append({"epoch": epoch, "train_loss": float(np.mean(losses)), "validation_pcc": val_metric["pcc"], "validation_mse": val_metric["mse"], "validation_mae": val_metric["mae"], "is_best": improved})
        pd.DataFrame(history).to_csv(out_dir / "metrics_history.csv", index=False)
        print(f"Epoch [{epoch:03d}/{args.epochs}]" + (" [BEST]" if improved else "") + f" loss={np.mean(losses):.5f} validation_PCC={val_metric['pcc']:.4f}", flush=True)
        if stale >= args.patience: print(f"Early stopped after {args.patience} epochs without validation PCC improvement.", flush=True); break
    model.load_state_dict(best_state)
    test_metric, prediction, details = evaluate(model, x, y, local_index, regional_index, local_position, regional_position, test_index, proteins, device, args.batch_size, save_details=True)
    np.save(out_dir / "test_predictions.npy", prediction); np.savez_compressed(out_dir / "test_v7_details.npz", test_index=test_index, **details); pd.DataFrame(test_metric.pop("per_protein")).to_csv(out_dir / "per_protein_test_metrics.csv", index=False)
    result = {"model": "TriScaleCenterProtein_V7_AnchoredCrossAttention", "setting": "same-slice strict spot holdout; whole-slice transductive RNA context only", "selection_protocol": "train-side validation PCC only", "pca_protocol": "fixed basis and mean fitted on fit proteins only", "dataset": args.dataset_name, "seed": args.seed, "n_spots": int(len(x)), "n_genes": int(x.shape[1]), "n_proteins": int(y.shape[1]), "n_fit": int(len(fit_index)), "n_validation": int(len(validation_index)), "n_test": int(len(test_index)), "local_neighbors": int(local_index.shape[1]), "regional_outer_neighbors": int(regional_index.shape[1]), "coordinate_scale": coordinate_scale, "program_components": int(n_components), "pca_explained_variance": float(pca.explained_variance_ratio_.sum()), "n_parameters": int(sum(p.numel() for p in model.parameters())), "best_epoch": int(best_epoch), "best_validation_pcc": float(best_pcc), "test": test_metric}
    (out_dir / "results.json").write_text(json.dumps(result, indent=2), encoding="utf-8"); print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
