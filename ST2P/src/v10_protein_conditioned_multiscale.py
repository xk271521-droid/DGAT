#!/usr/bin/env python3
"""V10: protein-conditioned multiscale V9 for cross-slice RNA-to-protein prediction.

This is an independent entrypoint: protected V9 code/results are never edited.
V10 keeps V9's corrected spatial neighbourhoods, shared RNA encoder, relative
position encoding, robustness perturbations, and source-only model selection.
It adds (i) protein-specific low-rank scale routing, (ii) protein-conditioned
neighbour cross-attention with shared keys/values, and (iii) a stable batchwise
protein-PCC loss alongside V9's MSE and auxiliary losses.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import shutil
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from sklearn.decomposition import PCA
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


def digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


class ProteinSpecificMLP(nn.Module):
    """29 independent nonlinear heads evaluated on protein-specific embeddings.

    Input is [batch, proteins, hidden], unlike V9's shared embedding.  The
    weights are still vectorised, so this does not instantiate 29 modules.
    """
    def __init__(self, hidden: int, proteins: int, width: int = 128):
        super().__init__()
        self.weight1 = nn.Parameter(torch.empty(proteins, hidden, width))
        self.bias1 = nn.Parameter(torch.empty(proteins, width))
        self.weight2 = nn.Parameter(torch.empty(proteins, width))
        self.bias2 = nn.Parameter(torch.empty(proteins))
        nn.init.uniform_(self.weight1, -1 / math.sqrt(hidden), 1 / math.sqrt(hidden))
        nn.init.uniform_(self.bias1, -1 / math.sqrt(hidden), 1 / math.sqrt(hidden))
        nn.init.uniform_(self.weight2, -1 / math.sqrt(width), 1 / math.sqrt(width))
        nn.init.uniform_(self.bias2, -1 / math.sqrt(width), 1 / math.sqrt(width))

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        # [B,P,H] x [P,H,W] -> [B,P,W] -> [B,P]
        h = F.leaky_relu(torch.einsum("bph,phw->bpw", z, self.weight1) + self.bias1)
        return (h * self.weight2).sum(dim=-1) + self.bias2


class ProteinAwareNeighbourPool(nn.Module):
    """Shared K/V attention, with one protein-conditioned query per protein.

    Query tokens have shape [B,P,H], whereas RNA neighbour keys/values remain
    [B,T,H].  MultiheadAttention shares all projections across proteins; only
    the low-dimensional protein embeddings make their queries distinct.
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

    def forward(self, query, context, key_padding_mask=None, need_attention=False):
        message, attention = self.attention(
            self.query_norm(query), self.context_norm(context), self.context_norm(context),
            key_padding_mask=key_padding_mask, need_weights=need_attention,
            average_attn_weights=False,
        )
        output = self.norm_attention(message)
        output = self.norm_ffn(output + self.dropout(self.ffn(output)))
        return output, attention


class ProteinConditionedV10(nn.Module):
    """V9 encoder plus protein-specific routing and neighbour attention."""
    def __init__(self, n_genes, n_proteins, hidden=256, heads=4, ffn_dim=256,
                 dropout=.3, protein_embedding_dim=32, router_rank=32):
        super().__init__()
        if hidden % heads:
            raise ValueError("hidden dimension must be divisible by attention heads")
        self.n_proteins = n_proteins
        self.rna_encoder = nn.Sequential(nn.Linear(n_genes, hidden), nn.ReLU(), nn.Dropout(dropout))
        self.center_norm = nn.LayerNorm(hidden)
        self.position_encoder = nn.Sequential(
            nn.Linear(3, hidden), nn.ReLU(), nn.Linear(hidden, hidden, bias=False),
        )
        self.protein_embedding = nn.Embedding(n_proteins, protein_embedding_dim)
        self.protein_query_projection = nn.Linear(protein_embedding_dim, hidden, bias=False)
        self.local_pool = ProteinAwareNeighbourPool(hidden, heads, ffn_dim, dropout)
        self.regional_pool = ProteinAwareNeighbourPool(hidden, heads, ffn_dim, dropout)

        # V9-style global scale preference plus a compact protein-specific
        # residual.  Parameter count is O(P * rank * 3), not O(P * H * 3).
        self.global_router = nn.Linear(hidden * 3, 3)
        self.router_context = nn.Linear(hidden * 3, router_rank)
        self.protein_router = nn.Parameter(torch.empty(n_proteins, router_rank, 3))
        self.protein_router_bias = nn.Parameter(torch.zeros(n_proteins, 3))
        nn.init.normal_(self.protein_router, mean=0.0, std=0.02)
        self.protein_head = ProteinSpecificMLP(hidden, n_proteins)

    def encode_neighbours(self, x, position):
        batch, tokens, genes = x.shape
        rna = self.rna_encoder(x.reshape(batch * tokens, genes)).reshape(batch, tokens, -1)
        return rna + self.position_encoder(position)

    def forward(self, center, local, regional, local_position, regional_position,
                local_key_padding_mask=None, regional_key_padding_mask=None, need_attention=False):
        batch = center.shape[0]
        z_center = self.center_norm(self.rna_encoder(center))                 # [B,H]
        protein_id = torch.arange(self.n_proteins, device=center.device)
        protein_query = z_center[:, None, :] + self.protein_query_projection(
            self.protein_embedding(protein_id)
        )[None, :, :]                                                          # [B,P,H]
        local_tokens = self.encode_neighbours(local, local_position)           # [B,4,H]
        regional_tokens = self.encode_neighbours(regional, regional_position)  # [B,4,H]
        z_local, local_attention = self.local_pool(
            protein_query, local_tokens, local_key_padding_mask, need_attention)
        z_regional, regional_attention = self.regional_pool(
            protein_query, regional_tokens, regional_key_padding_mask, need_attention)
        z_center_p = z_center[:, None, :].expand(batch, self.n_proteins, -1)
        levels = torch.stack((z_center_p, z_local, z_regional), dim=2)         # [B,P,3,H]
        router_input = torch.cat((z_center, z_local.mean(1), z_regional.mean(1)), dim=1)
        global_logits = self.global_router(router_input)[:, None, :]           # [B,1,3]
        context = self.router_context(router_input)                             # [B,R]
        protein_logits = torch.einsum("br,prs->bps", context, self.protein_router)
        weights = torch.softmax(global_logits + protein_logits + self.protein_router_bias[None], dim=-1)
        z_final = (levels * weights.unsqueeze(-1)).sum(dim=2)                  # [B,P,H]
        return {
            "final": self.protein_head(z_final),
            "center": self.protein_head(z_center_p),
            "local": self.protein_head(z_local),
            "regional": self.protein_head(z_regional),
            "protein_scale_attention": weights,
            "local_attention": local_attention,
            "regional_attention": regional_attention,
        }


def protein_pcc_loss(prediction, target, eps=1e-6):
    """One minus mean protein-wise Pearson over a minibatch, with safeguards."""
    pred = prediction - prediction.mean(dim=0, keepdim=True)
    truth = target - target.mean(dim=0, keepdim=True)
    denominator = torch.sqrt(pred.square().sum(0) * truth.square().sum(0)).clamp_min(eps)
    correlation = (pred * truth).sum(0) / denominator
    # The denominator clamp keeps gradients finite for near-constant proteins.
    valid = torch.isfinite(correlation) & (truth.square().sum(0) > eps)
    if not valid.any():
        return prediction.new_zeros(()), 0
    return 1.0 - correlation[valid].mean(), int(valid.sum().item())


def evaluate(model, x, y, local, regional, lp, rp, index, proteins, device, batch_size, details=False):
    from triscale_center_protein_v2_relpos import batch_inputs
    from triscale_center_protein_v9_dgat_compatible import dgat_metrics
    model.eval(); preds = []; scale = []; local_a = []; regional_a = []
    with torch.no_grad():
        for start in range(0, len(index), batch_size):
            idx = index[start:start + batch_size]
            out = model(*batch_inputs(x, local, regional, lp, rp, idx, device), need_attention=details)
            preds.append(out["final"].cpu().numpy())
            if details:
                scale.append(out["protein_scale_attention"].cpu().numpy())
                local_a.append(out["local_attention"].cpu().numpy())
                regional_a.append(out["regional_attention"].cpu().numpy())
    pred = np.concatenate(preds)
    diagnostic = None if not details else {
        "protein_scale_attention": np.concatenate(scale),
        "local_attention": np.concatenate(local_a),
        "regional_attention": np.concatenate(regional_a),
    }
    return dgat_metrics(pred, y[index], proteins), pred, diagnostic


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared_root", required=True); parser.add_argument("--out_dir", required=True)
    parser.add_argument("--experiment_dir", required=True); parser.add_argument("--target", default="A1_LN")
    parser.add_argument("--source_sections", nargs="+", required=True)
    parser.add_argument("--seed", type=int, default=42); parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--batch_size", type=int, default=64); parser.add_argument("--validation_ratio", type=float, default=.1)
    parser.add_argument("--hidden_dim", type=int, default=256); parser.add_argument("--heads", type=int, default=4)
    parser.add_argument("--ffn_dim", type=int, default=256); parser.add_argument("--dropout", type=float, default=.3)
    parser.add_argument("--protein_embedding_dim", type=int, default=32); parser.add_argument("--router_rank", type=int, default=32)
    parser.add_argument("--lr", type=float, default=5e-4); parser.add_argument("--weight_decay", type=float, default=1e-4)
    parser.add_argument("--aux_weight", type=float, default=.25); parser.add_argument("--pcc_weight", type=float, required=True)
    parser.add_argument("--gene_dropout_rate", type=float, default=.5); parser.add_argument("--neighbor_mask_rate", type=float, default=.25)
    parser.add_argument("--local_neighbors", type=int, default=4); parser.add_argument("--regional_neighbors", type=int, default=8)
    parser.add_argument("--program_components", type=int, default=9); parser.add_argument("--save_epoch_checkpoints", action="store_true")
    args = parser.parse_args()
    if not 0 <= args.pcc_weight or not 0 <= args.gene_dropout_rate < 1 or not 0 <= args.neighbor_mask_rate < 1:
        raise ValueError("Invalid loss/regularization weight")
    random.seed(args.seed); np.random.seed(args.seed); torch.manual_seed(args.seed)
    if torch.cuda.is_available(): torch.cuda.manual_seed_all(args.seed)
    torch.backends.cudnn.deterministic = True; torch.backends.cudnn.benchmark = False
    out = Path(args.out_dir)
    if out.exists(): raise FileExistsError(f"Refusing to overwrite {out}")
    out.mkdir(parents=True); (out / "config.json").write_text(json.dumps(vars(args), indent=2))
    if args.save_epoch_checkpoints: (out / "epoch_checkpoints").mkdir()
    shutil.copy2(__file__, out / "script_snapshot.py")
    sys.path.insert(0, str(Path(args.experiment_dir).resolve()))
    from v9_robust_source_subset_test_oracle import corrected_knn, dense, feature_dropout, panel, sample_key_padding_mask
    from triscale_center_protein_v2_relpos import relative_position, batch_inputs

    root = Path(args.prepared_root); manifest = json.loads((root / "prepared_manifest.json").read_text())
    if manifest.get("adt_clr_audit", {}).get("axis") != 0: raise ValueError("Require frozen Muon CLR axis=0 registry")
    names = [entry["dataset"] for entry in manifest["datasets"]]
    sources = list(args.source_sections)
    invalid = set(sources) - (set(names) - {args.target})
    if invalid or len(set(sources)) != len(sources): raise ValueError(f"Invalid/duplicate source sections: {sorted(invalid)}")
    genes, proteins = panel(root, names, "rna_dgat.h5ad"), panel(root, names, "adt_dgat_clr.h5ad")
    blocks, audit = {}, {}
    for name in [*sources, args.target]:
        rna = ad.read_h5ad(root / name / "rna_dgat.h5ad")[:, genes].copy()
        coords = np.asarray(rna.obsm["spatial"], dtype=np.float32)
        neighbours = corrected_knn(coords, args.regional_neighbors)
        local, regional = neighbours[:, :args.local_neighbors], neighbours[:, args.local_neighbors:]
        nearest = np.linalg.norm(coords[local[:, 0]] - coords, axis=1)
        scale = float(np.median(nearest[nearest > 0]))
        block = {"x": dense(rna.X), "local": local, "regional": regional,
                 "lp": relative_position(coords, local, scale), "rp": relative_position(coords, regional, scale),
                 "obs": np.asarray(rna.obs_names)}
        adt = ad.read_h5ad(root / name / "adt_dgat_clr.h5ad")[:, proteins].copy()
        if not rna.obs_names.equals(adt.obs_names): raise ValueError(f"RNA/ADT barcode mismatch: {name}")
        block["y"] = dense(adt.X); blocks[name] = block
        audit[name] = {"spots": len(rna), "neighbour_sha256": digest(neighbours), "coordinate_scale": scale}
    target = blocks[args.target]; target_index = np.arange(len(target["x"]))
    parts = [[] for _ in range(8)]; offset = 0; splits = {}
    for position, name in enumerate(sources):
        b = blocks[name]; order = np.arange(len(b["x"])); np.random.default_rng(args.seed + position).shuffle(order)
        n_val = max(1, int(round(len(order) * args.validation_ratio))); val, fit = np.sort(order[:n_val]), np.sort(order[n_val:])
        for dest, source in zip(parts[:6], [b["x"], b["y"], b["local"] + offset, b["regional"] + offset, b["lp"], b["rp"]]): dest.append(source)
        parts[6].append(fit + offset); parts[7].append(val + offset); splits[name] = {"fit": len(fit), "validation": len(val)}; offset += len(order)
    x, y, local, regional, lp, rp, fit, validation = map(np.concatenate, parts)
    np.savez_compressed(out / "source_split_indices.npz", fit=fit, validation=validation)
    pca = PCA(n_components=args.program_components, svd_solver="full", random_state=args.seed).fit(y[fit])
    np.save(out / "pca_basis_source_fit.npy", pca.components_); np.save(out / "pca_mean_source_fit.npy", pca.mean_)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = ProteinConditionedV10(len(genes), len(proteins), args.hidden_dim, args.heads, args.ffn_dim,
                                  args.dropout, args.protein_embedding_dim, args.router_rank).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    loader = DataLoader(TensorDataset(torch.from_numpy(fit)), batch_size=args.batch_size, shuffle=True,
                        generator=torch.Generator().manual_seed(args.seed))
    best_val, best_epoch, best_state = -np.inf, 0, None; best_test, best_test_epoch, best_test_state = -np.inf, 0, None; history = []
    print(f"V10 parameters={sum(v.numel() for v in model.parameters())}; sources={sources}; pcc_weight={args.pcc_weight}", flush=True)
    for epoch in range(1, args.epochs + 1):
        model.train(); total = []; final_losses = []; aux_losses = []; corr_losses = []; valid_counts = []
        for (idx_tensor,) in loader:
            idx = idx_tensor.numpy(); truth = torch.from_numpy(y[idx]).to(device)
            center, local_rna, regional_rna, local_pos, regional_pos = batch_inputs(x, local, regional, lp, rp, idx, device)
            output = model(feature_dropout(center, args.gene_dropout_rate), feature_dropout(local_rna, args.gene_dropout_rate),
                           feature_dropout(regional_rna, args.gene_dropout_rate), local_pos, regional_pos,
                           sample_key_padding_mask(len(idx), local_rna.shape[1], args.neighbor_mask_rate, device),
                           sample_key_padding_mask(len(idx), regional_rna.shape[1], args.neighbor_mask_rate, device))
            final_mse = F.mse_loss(output["final"], truth)
            auxiliary = sum(F.mse_loss(output[key], truth) for key in ("center", "local", "regional")) / 3
            corr, valid = protein_pcc_loss(output["final"], truth)
            loss = final_mse + args.aux_weight * auxiliary + args.pcc_weight * corr
            optimizer.zero_grad(set_to_none=True); loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step()
            total.append(loss.item()); final_losses.append(final_mse.item()); aux_losses.append(auxiliary.item()); corr_losses.append(corr.item()); valid_counts.append(valid)
        val, _, _ = evaluate(model, x, y, local, regional, lp, rp, validation, proteins, device, args.batch_size)
        test, _, _ = evaluate(model, target["x"], target["y"], target["local"], target["regional"], target["lp"], target["rp"], target_index, proteins, device, args.batch_size)
        if args.save_epoch_checkpoints: torch.save(model.state_dict(), out / "epoch_checkpoints" / f"epoch_{epoch:02d}.pth")
        source_best = val["pcc_pro_mean"] > best_val + 1e-6
        if source_best:
            best_val, best_epoch = val["pcc_pro_mean"], epoch; best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}; torch.save(best_state, out / "best_source_validation.pth")
        test_best = test["pcc_pro_mean"] > best_test + 1e-6
        if test_best:
            best_test, best_test_epoch = test["pcc_pro_mean"], epoch; best_test_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}; torch.save(best_test_state, out / "best_test_oracle.pth")
        row = {"epoch": epoch, "train_loss": float(np.mean(total)), "train_final_mse": float(np.mean(final_losses)),
               "train_aux_mse": float(np.mean(aux_losses)), "train_pcc_loss": float(np.mean(corr_losses)), "valid_proteins_per_batch": float(np.mean(valid_counts)),
               "source_validation_pcc_pro_mean": val["pcc_pro_mean"], "A1_test_pcc_pro_mean": test["pcc_pro_mean"], "A1_test_rmse_global": test["rmse_global"],
               "is_source_best": source_best, "is_test_oracle_best": test_best}
        history.append(row); pd.DataFrame(history).to_csv(out / "metrics_history.csv", index=False)
        print(f"Epoch {epoch:02d}/{args.epochs} loss={row['train_loss']:.4f} val_PCC={val['pcc_pro_mean']:.4f} A1_PCC={test['pcc_pro_mean']:.4f}" + (" [SOURCE_BEST]" if source_best else "") + (" [TEST_ORACLE]" if test_best else ""), flush=True)
    model.load_state_dict(best_state); strict, strict_pred, strict_details = evaluate(model, target["x"], target["y"], target["local"], target["regional"], target["lp"], target["rp"], target_index, proteins, device, args.batch_size, details=True)
    model.load_state_dict(best_test_state); oracle, oracle_pred, oracle_details = evaluate(model, target["x"], target["y"], target["local"], target["regional"], target["lp"], target["rp"], target_index, proteins, device, args.batch_size, details=True)
    for prefix, metric, prediction, detail in [("source_selected", strict, strict_pred, strict_details), ("test_oracle", oracle, oracle_pred, oracle_details)]:
        pd.DataFrame(metric.pop("per_protein")).to_csv(out / f"a1_{prefix}_per_protein_metrics.csv", index=False); np.save(out / f"a1_{prefix}_predictions.npy", prediction); np.savez_compressed(out / f"a1_{prefix}_attention.npz", **detail)
    result = {"model": "V10_ProteinConditionedMultiscale", "target": args.target, "sources": sources, "seed": args.seed, "epochs": args.epochs,
              "architecture": {"protein_specific_router": "global 3-scale logits plus low-rank protein residual", "protein_aware_neighbour_attention": "shared K/V MHA; protein embedding conditioned queries", "protein_specific_head": "vectorized 29 independent nonlinear heads", "protein_embedding_dim": args.protein_embedding_dim, "router_rank": args.router_rank},
              "loss": {"final_mse": 1.0, "protein_pcc": args.pcc_weight, "auxiliary_mse": args.aux_weight},
              "selection": "strict: D1 source-validation mean PCC-Pro; A1 evaluated each epoch only as diagnostic", "preprocessing": "frozen DGAT registry; Muon CLR axis=0", "common_panel": {"genes": len(genes), "proteins": proteins}, "source_split": splits, "audit": audit,
              "n_parameters": sum(v.numel() for v in model.parameters()), "best_epoch": best_epoch, "best_source_validation_pcc_pro_mean": best_val, "source_selected_test_A1": strict,
              "best_test_oracle_epoch": best_test_epoch, "best_test_oracle_pcc_pro_mean": best_test, "test_oracle_A1": oracle, "wall_seconds": time.monotonic()}
    (out / "results.json").write_text(json.dumps(result, indent=2)); print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
