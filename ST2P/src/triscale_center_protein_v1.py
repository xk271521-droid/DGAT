#!/usr/bin/env python3
"""TriScale-CenterProtein V1: a minimal same-slice spatial RNA-to-protein test.

Each scale is supervised to predict the *central* spot's protein panel:
center RNA -> center protein; local-neighbor RNA -> center protein; and
regional-neighbor RNA -> center protein.  Local and regional tokens are
pooled by center-conditioned attention, then the three level states undergo
one softmax fusion before the final protein decoder.

This is intentionally separate from V1-B/NPR-Net and does not use PCA,
protein graphs, prototypes, or target/test proteins during training.
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
from scipy.stats import pearsonr
from sklearn.neighbors import NearestNeighbors
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


def dense(x):
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x, dtype=np.float32)


def row_minmax(x: np.ndarray, low: float = 1e-8, high: float = 1.0) -> np.ndarray:
    lo, hi = x.min(axis=1, keepdims=True), x.max(axis=1, keepdims=True)
    out = low + (x - lo) / np.where(hi == lo, 1.0, hi - lo) * (high - low)
    out[(hi - lo).ravel() == 0] = low
    return np.where(np.isfinite(out), out, low).astype(np.float32)


def compute_metrics(pred: np.ndarray, target: np.ndarray, proteins: list[str]) -> dict:
    per_protein = []
    for j, protein in enumerate(proteins):
        a, b = pred[:, j], target[:, j]
        pcc = np.nan if np.std(a) == 0 or np.std(b) == 0 else float(pearsonr(a, b)[0])
        per_protein.append({
            "protein": protein,
            "pcc": pcc,
            "mse": float(np.mean((a - b) ** 2)),
            "mae": float(np.mean(np.abs(a - b))),
        })
    return {
        "pcc": float(np.nanmean([r["pcc"] for r in per_protein])),
        "pcc_std": float(np.nanstd([r["pcc"] for r in per_protein])),
        "mse": float(np.mean((pred - target) ** 2)),
        "mae": float(np.mean(np.abs(pred - target))),
        "per_protein": per_protein,
    }


class TriScaleCenterProtein(nn.Module):
    """One shared RNA encoder, within-level attention, and one level fusion."""

    def __init__(self, n_genes: int, n_proteins: int, hidden: int, heads: int, dropout: float):
        super().__init__()
        if hidden % heads:
            raise ValueError("hidden dimension must be divisible by attention heads")
        self.rna_encoder = nn.Sequential(
            nn.Linear(n_genes, hidden), nn.ReLU(), nn.Dropout(dropout)
        )
        self.local_pool = nn.MultiheadAttention(hidden, heads, dropout=dropout, batch_first=True)
        self.regional_pool = nn.MultiheadAttention(hidden, heads, dropout=dropout, batch_first=True)
        self.local_norm = nn.LayerNorm(hidden)
        self.regional_norm = nn.LayerNorm(hidden)
        self.level_gate = nn.Sequential(nn.Linear(hidden * 3, 3), nn.Softmax(dim=1))
        # A shared decoder keeps the three auxiliary tasks in a common protein space.
        self.protein_decoder = nn.Linear(hidden, n_proteins)

    def encode_tokens(self, x: torch.Tensor) -> torch.Tensor:
        """Encode [B, L, G] RNA tokens with the shared RNA encoder."""
        batch, tokens, genes = x.shape
        return self.rna_encoder(x.reshape(batch * tokens, genes)).reshape(batch, tokens, -1)

    def forward(self, center: torch.Tensor, local: torch.Tensor, regional: torch.Tensor, need_attention: bool = False):
        # center [B,G]; local [B,6,G]; regional [B,18,G]
        h_center = self.rna_encoder(center)  # [B,H]
        q = h_center.unsqueeze(1)            # [B,1,H]
        h_local = self.encode_tokens(local)  # [B,6,H]
        h_regional = self.encode_tokens(regional)  # [B,18,H]

        # The center queries each level but the returned state is formed from
        # neighbor values only; a local/regional branch cannot copy center RNA.
        z_local, a_local = self.local_pool(q, h_local, h_local, need_weights=need_attention,
                                            average_attn_weights=False)
        z_regional, a_regional = self.regional_pool(q, h_regional, h_regional, need_weights=need_attention,
                                                      average_attn_weights=False)
        z_local = self.local_norm(z_local.squeeze(1))
        z_regional = self.regional_norm(z_regional.squeeze(1))

        levels = torch.stack((h_center, z_local, z_regional), dim=1)  # [B,3,H]
        alpha = self.level_gate(torch.cat((h_center, z_local, z_regional), dim=1))  # [B,3]
        z_fused = (levels * alpha.unsqueeze(-1)).sum(dim=1)  # [B,H]

        pred_center = self.protein_decoder(h_center)
        pred_local = self.protein_decoder(z_local)
        pred_regional = self.protein_decoder(z_regional)
        pred_final = self.protein_decoder(z_fused)
        return {
            "final": pred_final,
            "center": pred_center,
            "local": pred_local,
            "regional": pred_regional,
            "level_weights": alpha,
            "local_attention": a_local,
            "regional_attention": a_regional,
        }


def batch_inputs(x, local_index, regional_index, index, device):
    index = np.asarray(index, dtype=np.int64)
    center = torch.from_numpy(x[index]).to(device, non_blocking=True)
    local = torch.from_numpy(x[local_index[index]]).to(device, non_blocking=True)
    regional = torch.from_numpy(x[regional_index[index]]).to(device, non_blocking=True)
    return center, local, regional


def evaluate(model, x, y, local_index, regional_index, index, proteins, device, batch_size, save_details=False):
    model.eval()
    prediction, levels, local_attn, regional_attn = [], [], [], []
    with torch.no_grad():
        for start in range(0, len(index), batch_size):
            batch_index = index[start:start + batch_size]
            center, local, regional = batch_inputs(x, local_index, regional_index, batch_index, device)
            out = model(center, local, regional, need_attention=save_details)
            prediction.append(out["final"].cpu().numpy())
            if save_details:
                levels.append(out["level_weights"].cpu().numpy())
                local_attn.append(out["local_attention"].cpu().numpy())
                regional_attn.append(out["regional_attention"].cpu().numpy())
    pred = np.concatenate(prediction, axis=0)
    details = None
    if save_details:
        details = {
            "level_weights": np.concatenate(levels, axis=0),
            "local_attention": np.concatenate(local_attn, axis=0),
            "regional_attention": np.concatenate(regional_attn, axis=0),
        }
    return compute_metrics(pred, y[index], proteins), pred, details


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna_h5ad", required=True)
    parser.add_argument("--adt_h5ad", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--dataset_name", default="Human Lymph Node A1")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--train_ratio", type=float, default=0.9)
    parser.add_argument("--validation_ratio", type=float, default=0.1)
    parser.add_argument("--local_neighbors", type=int, default=6)
    parser.add_argument("--regional_neighbors", type=int, default=24)
    parser.add_argument("--epochs", type=int, default=500)
    parser.add_argument("--patience", type=int, default=50)
    parser.add_argument("--batch_size", type=int, default=128)
    parser.add_argument("--hidden_dim", type=int, default=256)
    parser.add_argument("--heads", type=int, default=4)
    parser.add_argument("--dropout", type=float, default=0.3)
    parser.add_argument("--lr", type=float, default=5e-4)
    parser.add_argument("--weight_decay", type=float, default=1e-4)
    parser.add_argument("--aux_weight", type=float, default=0.25)
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    out_dir = Path(args.out_dir)
    if out_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing result directory: {out_dir}")
    out_dir.mkdir(parents=True)
    (out_dir / "config.json").write_text(json.dumps(vars(args), indent=2), encoding="utf-8")

    rna, adt = ad.read_h5ad(args.rna_h5ad), ad.read_h5ad(args.adt_h5ad)
    shared = rna.obs_names.intersection(adt.obs_names, sort=False)
    rna, adt = rna[shared].copy(), adt[shared].copy()
    if "spatial" not in rna.obsm:
        raise KeyError("RNA h5ad must contain obsm['spatial']")
    x, y = row_minmax(dense(rna.X)), row_minmax(dense(adt.X))
    keep = np.var(y, axis=0) > 0
    y, proteins = y[:, keep], [str(p) for p in adt.var_names[keep]]
    coords = np.asarray(rna.obsm["spatial"], dtype=np.float32)

    neighbor_model = NearestNeighbors(n_neighbors=args.regional_neighbors + 1).fit(coords)
    all_neighbors = neighbor_model.kneighbors(return_distance=False)[:, 1:]
    local_index = all_neighbors[:, :args.local_neighbors]
    # Use a non-overlapping outer ring so regional input does not duplicate local neighbors.
    regional_index = all_neighbors[:, args.local_neighbors:args.regional_neighbors]
    if regional_index.shape[1] == 0:
        raise ValueError("regional_neighbors must exceed local_neighbors")

    rng = np.random.default_rng(args.seed)
    shuffled = np.arange(len(x)); rng.shuffle(shuffled)
    outer_cut = int(len(shuffled) * args.train_ratio)
    outer_train, test_index = np.sort(shuffled[:outer_cut]), np.sort(shuffled[outer_cut:])
    inner = outer_train.copy(); rng = np.random.default_rng(args.seed); rng.shuffle(inner)
    val_n = max(1, int(round(len(inner) * args.validation_ratio)))
    validation_index, fit_index = np.sort(inner[:val_n]), np.sort(inner[val_n:])
    np.savez_compressed(out_dir / "split_indices.npz", fit=fit_index, validation=validation_index, test=test_index)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = TriScaleCenterProtein(x.shape[1], y.shape[1], args.hidden_dim, args.heads, args.dropout).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    loader = DataLoader(TensorDataset(torch.from_numpy(fit_index)), batch_size=args.batch_size,
                        shuffle=True, generator=torch.Generator().manual_seed(args.seed))
    print(f"TriScaleCenterProtein parameters: {sum(p.numel() for p in model.parameters()):,}", flush=True)
    print(f"Spots fit/validation/test: {len(fit_index)}/{len(validation_index)}/{len(test_index)}; "
          f"Local/outer-regional tokens: {local_index.shape[1]}/{regional_index.shape[1]}", flush=True)

    best_pcc, best_epoch, stale, best_state, history = -np.inf, 0, 0, None, []
    mse = nn.MSELoss()
    for epoch in range(1, args.epochs + 1):
        model.train(); epoch_losses = []
        for (batch_index,) in loader:
            idx = batch_index.numpy()
            center, local, regional = batch_inputs(x, local_index, regional_index, idx, device)
            target = torch.from_numpy(y[idx]).to(device, non_blocking=True)
            out = model(center, local, regional)
            auxiliary = (mse(out["center"], target) + mse(out["local"], target) + mse(out["regional"], target)) / 3.0
            loss = mse(out["final"], target) + args.aux_weight * auxiliary
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            epoch_losses.append(float(loss.detach().cpu()))

        val_metric, _, _ = evaluate(model, x, y, local_index, regional_index, validation_index, proteins,
                                     device, args.batch_size, save_details=False)
        improved = val_metric["pcc"] > best_pcc + 1e-6
        if improved:
            best_pcc, best_epoch, stale = val_metric["pcc"], epoch, 0
            best_state = deepcopy({k: v.detach().cpu() for k, v in model.state_dict().items()})
            torch.save(best_state, out_dir / "best_validation.pth")
        else:
            stale += 1
        history.append({"epoch": epoch, "train_loss": float(np.mean(epoch_losses)),
                        "validation_pcc": val_metric["pcc"], "validation_mse": val_metric["mse"],
                        "validation_mae": val_metric["mae"], "is_best": improved})
        pd.DataFrame(history).to_csv(out_dir / "metrics_history.csv", index=False)
        print(f"Epoch [{epoch:03d}/{args.epochs}]" + (" [BEST]" if improved else "") +
              f" loss={np.mean(epoch_losses):.5f} validation_PCC={val_metric['pcc']:.4f}", flush=True)
        if stale >= args.patience:
            print(f"Early stopped after {args.patience} epochs without validation PCC improvement.", flush=True)
            break

    model.load_state_dict(best_state)
    test_metric, pred, details = evaluate(model, x, y, local_index, regional_index, test_index, proteins,
                                          device, args.batch_size, save_details=True)
    np.save(out_dir / "test_predictions.npy", pred)
    np.savez_compressed(out_dir / "test_attention.npz", test_index=test_index, **details)
    pd.DataFrame(test_metric.pop("per_protein")).to_csv(out_dir / "per_protein_test_metrics.csv", index=False)
    result = {
        "model": "TriScaleCenterProtein_V1",
        "setting": "same-slice strict spot holdout; whole-slice transductive RNA context only",
        "selection_protocol": "train-side validation PCC only",
        "dataset": args.dataset_name,
        "seed": args.seed,
        "n_spots": int(len(x)), "n_genes": int(x.shape[1]), "n_proteins": int(y.shape[1]),
        "n_fit": int(len(fit_index)), "n_validation": int(len(validation_index)), "n_test": int(len(test_index)),
        "local_neighbors": int(local_index.shape[1]), "regional_outer_neighbors": int(regional_index.shape[1]),
        "n_parameters": int(sum(p.numel() for p in model.parameters())), "best_epoch": int(best_epoch),
        "best_validation_pcc": float(best_pcc), "test": test_metric,
    }
    (out_dir / "results.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
