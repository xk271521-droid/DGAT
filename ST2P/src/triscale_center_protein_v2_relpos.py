#!/usr/bin/env python3
"""TriScale-CenterProtein V2: V1 plus center-relative spatial position tokens.

The only architectural change from V1 is that every local/regional RNA token
receives an embedding of (delta_x, delta_y, Euclidean_distance) relative to
the center spot before center-guided attention pooling.
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
from sklearn.neighbors import NearestNeighbors
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

from triscale_center_protein_v1 import compute_metrics, dense, row_minmax


class TriScaleCenterProteinRelPos(nn.Module):
    """V1 model with learned relative-position embeddings for neighbor tokens.

    ``center_layernorm=False`` exactly retains the completed V2 architecture.
    ``relation_pool=False`` retains the completed V2/V3 multi-head pooling path.
    The two optional switches are reserved for isolated later ablations.
    """

    def __init__(self, n_genes: int, n_proteins: int, hidden: int, heads: int, dropout: float,
                 center_layernorm: bool = False, relation_pool: bool = False,
                 relation_hidden: int = 64):
        super().__init__()
        if hidden % heads:
            raise ValueError("hidden dimension must be divisible by attention heads")
        self.rna_encoder = nn.Sequential(nn.Linear(n_genes, hidden), nn.ReLU(), nn.Dropout(dropout))
        self.center_norm = nn.LayerNorm(hidden) if center_layernorm else nn.Identity()
        self.relation_pool = relation_pool
        self.position_encoder = nn.Sequential(
            nn.Linear(3, hidden), nn.ReLU(), nn.Linear(hidden, hidden, bias=False)
        )
        if relation_pool:
            # A small additive scorer is deliberately used instead of another
            # transformer: it makes every neighbour's center/content/position
            # relation directly responsible for its aggregation weight.
            self.local_center_proj = nn.Linear(hidden, relation_hidden, bias=False)
            self.local_token_proj = nn.Linear(hidden, relation_hidden, bias=False)
            self.local_position_proj = nn.Linear(hidden, relation_hidden, bias=False)
            self.local_relation_score = nn.Linear(relation_hidden, 1, bias=False)
            self.regional_center_proj = nn.Linear(hidden, relation_hidden, bias=False)
            self.regional_token_proj = nn.Linear(hidden, relation_hidden, bias=False)
            self.regional_position_proj = nn.Linear(hidden, relation_hidden, bias=False)
            self.regional_relation_score = nn.Linear(relation_hidden, 1, bias=False)
        else:
            self.local_pool = nn.MultiheadAttention(hidden, heads, dropout=dropout, batch_first=True)
            self.regional_pool = nn.MultiheadAttention(hidden, heads, dropout=dropout, batch_first=True)
        self.local_norm = nn.LayerNorm(hidden)
        self.regional_norm = nn.LayerNorm(hidden)
        self.level_gate = nn.Sequential(nn.Linear(hidden * 3, 3), nn.Softmax(dim=1))
        self.protein_decoder = nn.Linear(hidden, n_proteins)

    def encode_tokens(self, x: torch.Tensor, relative_position: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        batch, tokens, genes = x.shape
        rna = self.rna_encoder(x.reshape(batch * tokens, genes)).reshape(batch, tokens, -1)
        return rna, self.position_encoder(relative_position)

    def relation_aggregate(self, center: torch.Tensor, tokens: torch.Tensor, positions: torch.Tensor,
                           center_proj: nn.Linear, token_proj: nn.Linear, position_proj: nn.Linear,
                           score_proj: nn.Linear) -> tuple[torch.Tensor, torch.Tensor]:
        relation = torch.tanh(center_proj(center).unsqueeze(1) + token_proj(tokens) + position_proj(positions))
        weight = torch.softmax(score_proj(relation).squeeze(-1), dim=1)
        return (weight.unsqueeze(-1) * (tokens + positions)).sum(dim=1), weight

    def forward(self, center, local, regional, local_position, regional_position, need_attention: bool = False):
        h_center = self.center_norm(self.rna_encoder(center))
        q = h_center.unsqueeze(1)
        local_rna, local_pos = self.encode_tokens(local, local_position)
        regional_rna, regional_pos = self.encode_tokens(regional, regional_position)
        if self.relation_pool:
            z_local, a_local = self.relation_aggregate(
                h_center, local_rna, local_pos, self.local_center_proj, self.local_token_proj,
                self.local_position_proj, self.local_relation_score)
            z_regional, a_regional = self.relation_aggregate(
                h_center, regional_rna, regional_pos, self.regional_center_proj, self.regional_token_proj,
                self.regional_position_proj, self.regional_relation_score)
        else:
            h_local, h_regional = local_rna + local_pos, regional_rna + regional_pos
            z_local, a_local = self.local_pool(q, h_local, h_local, need_weights=need_attention,
                                                average_attn_weights=False)
            z_regional, a_regional = self.regional_pool(q, h_regional, h_regional, need_weights=need_attention,
                                                          average_attn_weights=False)
            z_local, z_regional = z_local.squeeze(1), z_regional.squeeze(1)
        z_local = self.local_norm(z_local)
        z_regional = self.regional_norm(z_regional)
        levels = torch.stack((h_center, z_local, z_regional), dim=1)
        alpha = self.level_gate(torch.cat((h_center, z_local, z_regional), dim=1))
        z_fused = (levels * alpha.unsqueeze(-1)).sum(dim=1)
        return {
            "final": self.protein_decoder(z_fused),
            "center": self.protein_decoder(h_center),
            "local": self.protein_decoder(z_local),
            "regional": self.protein_decoder(z_regional),
            "level_weights": alpha,
            "local_attention": a_local,
            "regional_attention": a_regional,
        }


def relative_position(coords: np.ndarray, neighbors: np.ndarray, scale: float) -> np.ndarray:
    delta = (coords[neighbors] - coords[:, None, :]) / scale
    distance = np.linalg.norm(delta, axis=-1, keepdims=True)
    return np.concatenate((delta, distance), axis=-1).astype(np.float32)


def batch_inputs(x, local_index, regional_index, local_position, regional_position, index, device):
    index = np.asarray(index, dtype=np.int64)
    return (
        torch.from_numpy(x[index]).to(device, non_blocking=True),
        torch.from_numpy(x[local_index[index]]).to(device, non_blocking=True),
        torch.from_numpy(x[regional_index[index]]).to(device, non_blocking=True),
        torch.from_numpy(local_position[index]).to(device, non_blocking=True),
        torch.from_numpy(regional_position[index]).to(device, non_blocking=True),
    )


def evaluate(model, x, y, local_index, regional_index, local_position, regional_position, index, proteins, device, batch_size, save_details=False):
    model.eval(); prediction, levels, local_attn, regional_attn = [], [], [], []
    with torch.no_grad():
        for start in range(0, len(index), batch_size):
            idx = index[start:start + batch_size]
            out = model(*batch_inputs(x, local_index, regional_index, local_position, regional_position, idx, device),
                        need_attention=save_details)
            prediction.append(out["final"].cpu().numpy())
            if save_details:
                levels.append(out["level_weights"].cpu().numpy())
                local_attn.append(out["local_attention"].cpu().numpy())
                regional_attn.append(out["regional_attention"].cpu().numpy())
    details = None
    if save_details:
        details = {"level_weights": np.concatenate(levels), "local_attention": np.concatenate(local_attn),
                   "regional_attention": np.concatenate(regional_attn)}
    pred = np.concatenate(prediction)
    return compute_metrics(pred, y[index], proteins), pred, details


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna_h5ad", required=True); ap.add_argument("--adt_h5ad", required=True); ap.add_argument("--out_dir", required=True)
    ap.add_argument("--dataset_name", default="Human Lymph Node A1"); ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--train_ratio", type=float, default=0.9); ap.add_argument("--validation_ratio", type=float, default=0.1)
    ap.add_argument("--local_neighbors", type=int, default=6); ap.add_argument("--regional_neighbors", type=int, default=24)
    ap.add_argument("--epochs", type=int, default=500); ap.add_argument("--patience", type=int, default=50)
    ap.add_argument("--batch_size", type=int, default=64); ap.add_argument("--hidden_dim", type=int, default=256)
    ap.add_argument("--heads", type=int, default=4); ap.add_argument("--dropout", type=float, default=0.3)
    ap.add_argument("--lr", type=float, default=5e-4); ap.add_argument("--weight_decay", type=float, default=1e-4)
    ap.add_argument("--aux_weight", type=float, default=0.25)
    ap.add_argument("--center_layernorm", action="store_true",
                    help="V3 ablation: align the Center branch scale with Local and Regional branches.")
    ap.add_argument("--relation_pool", action="store_true",
                    help="V4: use center-conditioned RNA-position relation scoring for neighbour pooling.")
    ap.add_argument("--relation_hidden", type=int, default=64)
    args = ap.parse_args()
    random.seed(args.seed); np.random.seed(args.seed); torch.manual_seed(args.seed)
    if torch.cuda.is_available(): torch.cuda.manual_seed_all(args.seed)
    torch.backends.cudnn.deterministic = True; torch.backends.cudnn.benchmark = False

    out_dir = Path(args.out_dir)
    if out_dir.exists(): raise FileExistsError(f"Refusing to overwrite existing result directory: {out_dir}")
    out_dir.mkdir(parents=True)
    (out_dir / "config.json").write_text(json.dumps(vars(args), indent=2), encoding="utf-8")
    rna, adt = ad.read_h5ad(args.rna_h5ad), ad.read_h5ad(args.adt_h5ad)
    shared = rna.obs_names.intersection(adt.obs_names, sort=False); rna, adt = rna[shared].copy(), adt[shared].copy()
    if "spatial" not in rna.obsm: raise KeyError("RNA h5ad must contain obsm['spatial']")
    x, y = row_minmax(dense(rna.X)), row_minmax(dense(adt.X))
    keep = np.var(y, axis=0) > 0; y, proteins = y[:, keep], [str(p) for p in adt.var_names[keep]]
    coords = np.asarray(rna.obsm["spatial"], dtype=np.float32)
    all_neighbors = NearestNeighbors(n_neighbors=args.regional_neighbors + 1).fit(coords).kneighbors(return_distance=False)[:, 1:]
    local_index, regional_index = all_neighbors[:, :args.local_neighbors], all_neighbors[:, args.local_neighbors:args.regional_neighbors]
    if regional_index.shape[1] == 0: raise ValueError("regional_neighbors must exceed local_neighbors")
    nearest_distance = np.linalg.norm(coords[local_index[:, 0]] - coords, axis=1)
    coordinate_scale = float(np.median(nearest_distance[nearest_distance > 0]))
    if not np.isfinite(coordinate_scale) or coordinate_scale <= 0: raise ValueError("invalid coordinate scale")
    local_position = relative_position(coords, local_index, coordinate_scale)
    regional_position = relative_position(coords, regional_index, coordinate_scale)

    rng = np.random.default_rng(args.seed); shuffled = np.arange(len(x)); rng.shuffle(shuffled)
    outer_cut = int(len(shuffled) * args.train_ratio); outer_train, test_index = np.sort(shuffled[:outer_cut]), np.sort(shuffled[outer_cut:])
    inner = outer_train.copy(); rng = np.random.default_rng(args.seed); rng.shuffle(inner)
    val_n = max(1, int(round(len(inner) * args.validation_ratio))); validation_index, fit_index = np.sort(inner[:val_n]), np.sort(inner[val_n:])
    np.savez_compressed(out_dir / "split_indices.npz", fit=fit_index, validation=validation_index, test=test_index)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = TriScaleCenterProteinRelPos(
        x.shape[1], y.shape[1], args.hidden_dim, args.heads, args.dropout,
        center_layernorm=args.center_layernorm, relation_pool=args.relation_pool,
        relation_hidden=args.relation_hidden,
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    loader = DataLoader(TensorDataset(torch.from_numpy(fit_index)), batch_size=args.batch_size, shuffle=True,
                        generator=torch.Generator().manual_seed(args.seed))
    model_name = ("TriScaleCenterProtein_V4_RelationPool" if args.relation_pool else
                  "TriScaleCenterProtein_V3_CenterNorm" if args.center_layernorm else
                  "TriScaleCenterProtein_V2_RelPos")
    print(f"{model_name} parameters: {sum(p.numel() for p in model.parameters()):,}; coordinate_scale={coordinate_scale:.4f}", flush=True)
    mse = nn.MSELoss(); best_pcc, best_epoch, stale, best_state, history = -np.inf, 0, 0, None, []
    for epoch in range(1, args.epochs + 1):
        model.train(); losses = []
        for (batch_index,) in loader:
            idx = batch_index.numpy(); target = torch.from_numpy(y[idx]).to(device, non_blocking=True)
            out = model(*batch_inputs(x, local_index, regional_index, local_position, regional_position, idx, device))
            auxiliary = (mse(out["center"], target) + mse(out["local"], target) + mse(out["regional"], target)) / 3.0
            loss = mse(out["final"], target) + args.aux_weight * auxiliary
            optimizer.zero_grad(set_to_none=True); loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step()
            losses.append(float(loss.detach().cpu()))
        val_metric, _, _ = evaluate(model, x, y, local_index, regional_index, local_position, regional_position,
                                     validation_index, proteins, device, args.batch_size)
        improved = val_metric["pcc"] > best_pcc + 1e-6
        if improved:
            best_pcc, best_epoch, stale = val_metric["pcc"], epoch, 0
            best_state = deepcopy({k: v.detach().cpu() for k, v in model.state_dict().items()})
            torch.save(best_state, out_dir / "best_validation.pth")
        else: stale += 1
        history.append({"epoch": epoch, "train_loss": float(np.mean(losses)), "validation_pcc": val_metric["pcc"],
                        "validation_mse": val_metric["mse"], "validation_mae": val_metric["mae"], "is_best": improved})
        pd.DataFrame(history).to_csv(out_dir / "metrics_history.csv", index=False)
        print(f"Epoch [{epoch:03d}/{args.epochs}]" + (" [BEST]" if improved else "") +
              f" loss={np.mean(losses):.5f} validation_PCC={val_metric['pcc']:.4f}", flush=True)
        if stale >= args.patience:
            print(f"Early stopped after {args.patience} epochs without validation PCC improvement.", flush=True); break
    model.load_state_dict(best_state)
    test_metric, pred, details = evaluate(model, x, y, local_index, regional_index, local_position, regional_position,
                                          test_index, proteins, device, args.batch_size, save_details=True)
    np.save(out_dir / "test_predictions.npy", pred); np.savez_compressed(out_dir / "test_attention.npz", test_index=test_index, **details)
    pd.DataFrame(test_metric.pop("per_protein")).to_csv(out_dir / "per_protein_test_metrics.csv", index=False)
    result = {"model": model_name, "setting": "same-slice strict spot holdout; whole-slice transductive RNA context only",
              "selection_protocol": "train-side validation PCC only", "dataset": args.dataset_name, "seed": args.seed,
              "n_spots": int(len(x)), "n_genes": int(x.shape[1]), "n_proteins": int(y.shape[1]),
              "n_fit": int(len(fit_index)), "n_validation": int(len(validation_index)), "n_test": int(len(test_index)),
              "local_neighbors": int(local_index.shape[1]), "regional_outer_neighbors": int(regional_index.shape[1]),
              "coordinate_scale": coordinate_scale, "n_parameters": int(sum(p.numel() for p in model.parameters())),
              "best_epoch": int(best_epoch), "best_validation_pcc": float(best_pcc), "test": test_metric}
    (out_dir / "results.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
