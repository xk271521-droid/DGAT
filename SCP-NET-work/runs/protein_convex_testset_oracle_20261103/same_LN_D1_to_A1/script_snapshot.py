#!/usr/bin/env python3
"""SCN-Pro: Slice-Conditioned Niche-to-Protein Network.

This independent pilot preserves V10's protected local multi-scale niche
encoder, relative coordinates, protein-specific scale router, robust input
perturbations and source-only checkpoint selection.  It adds only RNA-visible
cross-slice mechanisms: (1) whole-slice RNA context, (2) per-slice input
canonicalization, and (3) a low-rank context adapter in the protein decoder.
Target ADT is never used by these mechanisms.
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



class CognateAnchorModule(nn.Module):
    """Vectorized biological prior anchor linking multi-scale cognate RNA to predictions."""
    def __init__(self, n_proteins: int, hidden: int):
        super().__init__()
        self.n_proteins = n_proteins
        # 3 multi-scale inputs: [center, local_mean, regional_mean]
        self.latent_proj = nn.Parameter(torch.zeros(n_proteins, 3, hidden))
        self.direct_weight = nn.Parameter(torch.zeros(n_proteins, 3))
        self.direct_bias = nn.Parameter(torch.zeros(n_proteins))
        # Learnable gating logit per protein initialized to -2.0 (sigmoid ~ 0.12)
        self.gate_logit = nn.Parameter(torch.full((n_proteins,), -2.0))
        # Per-protein learnable convex blend between dynamic attention and static prior
        self.alpha_logit = nn.Parameter(torch.full((n_proteins,), -2.0))
        with torch.no_grad():
            self.direct_weight[:, 0] = 0.5
            self.direct_weight[:, 1] = 0.3
            self.direct_weight[:, 2] = 0.2

    def forward(self, u: torch.Tensor, z: torch.Tensor, mask: torch.Tensor, scale_weights: torch.Tensor = None):
        # u: [B, P, 3], z: [B, P, H], mask: [P]
        u_masked = u * mask[None, :, None]
        z_inj = torch.einsum("bpc,pch->bph", u_masked, self.latent_proj)
        z_out = z + z_inj
        gate = torch.sigmoid(self.gate_logit) * mask
        # Protein-Specific Convex Routing:
        if scale_weights is not None:
            alpha = torch.sigmoid(self.alpha_logit)[None, :, None]
            sw = alpha * scale_weights + (1.0 - alpha) * self.direct_weight[None, :, :]
        else:
            sw = self.direct_weight[None, :, :]
        direct_pred = (u_masked * sw).sum(dim=-1) + self.direct_bias[None, :]
        return z_out, direct_pred, gate


class SliceConditionedNicheProtein(nn.Module):
    """V10 local niche prediction conditioned on RNA-only slice context."""
    def __init__(self, n_genes, n_proteins, hidden=256, heads=4, ffn_dim=256,
                 dropout=.3, protein_embedding_dim=32, router_rank=32,
                 context_input_dim=67, context_rank=8, protein_specific_router=True,
                 protein_aware_attention=True, use_slice_context=True,
                 use_canonicalization=True, use_context_adapter=True,
                 cognate_indices=None, has_cognate=None, use_cognate_anchor=True):
        super().__init__()
        if hidden % heads:
            raise ValueError("hidden dimension must be divisible by attention heads")
        self.n_proteins = n_proteins
        self.protein_specific_router_enabled = protein_specific_router
        self.protein_aware_attention_enabled = protein_aware_attention
        self.use_slice_context = use_slice_context
        self.use_canonicalization = use_canonicalization
        self.use_context_adapter = use_context_adapter
        self.use_cognate_anchor = use_cognate_anchor
        if cognate_indices is not None and has_cognate is not None:
            self.register_buffer("cognate_indices", torch.as_tensor(cognate_indices, dtype=torch.long))
            self.register_buffer("has_cognate", torch.as_tensor(has_cognate, dtype=torch.float32))
        else:
            self.register_buffer("cognate_indices", torch.zeros(n_proteins, dtype=torch.long))
            self.register_buffer("has_cognate", torch.zeros(n_proteins, dtype=torch.float32))
        self.rna_encoder = nn.Sequential(nn.Linear(n_genes, hidden), nn.ReLU(), nn.Dropout(dropout))
        self.center_norm = nn.LayerNorm(hidden)
        self.position_encoder = nn.Sequential(
            nn.Linear(3, hidden), nn.ReLU(), nn.Linear(hidden, hidden, bias=False),
        )
        if protein_aware_attention:
            self.protein_embedding = nn.Embedding(n_proteins, protein_embedding_dim)
            self.protein_query_projection = nn.Linear(protein_embedding_dim, hidden, bias=False)
        self.local_pool = ProteinAwareNeighbourPool(hidden, heads, ffn_dim, dropout)
        self.regional_pool = ProteinAwareNeighbourPool(hidden, heads, ffn_dim, dropout)

        # V9-style global scale preference plus a compact protein-specific
        # residual.  Parameter count is O(P * rank * 3), not O(P * H * 3).
        self.global_router = nn.Linear(hidden * 3, 3)
        if protein_specific_router:
            self.router_context = nn.Linear(hidden * 3, router_rank)
            self.protein_router = nn.Parameter(torch.empty(n_proteins, router_rank, 3))
            self.protein_router_bias = nn.Parameter(torch.zeros(n_proteins, 3))
            nn.init.normal_(self.protein_router, mean=0.0, std=0.02)
        # The context descriptor is constructed from source-fitted RNA PCA
        # statistics of the complete slice (and never from protein values).
        self.context_encoder = nn.Sequential(
            nn.Linear(context_input_dim, 64), nn.GELU(), nn.Dropout(dropout), nn.Linear(64, hidden), nn.Tanh(),
        )
        self.context_router = nn.Linear(hidden, hidden * 3, bias=False)
        nn.init.zeros_(self.context_router.weight)
        # Start from the protected V10 input path and learn only as much
        # slice-wise canonicalization as source data supports.
        self.canonical_logit = nn.Parameter(torch.tensor(-2.2))
        self.context_to_rank = nn.Linear(hidden, context_rank, bias=False)
        self.context_adapter_basis = nn.Parameter(torch.zeros(n_proteins, context_rank, hidden))
        self.protein_head = ProteinSpecificMLP(hidden, n_proteins)
        if self.use_cognate_anchor:
            self.cognate_anchor = CognateAnchorModule(n_proteins, hidden)

    def encode_neighbours(self, x, position):
        batch, tokens, genes = x.shape
        rna = self.rna_encoder(x.reshape(batch * tokens, genes)).reshape(batch, tokens, -1)
        return rna + self.position_encoder(position)

    def context_adapt(self, level, context):
        """Add a bounded protein-specific residual initialized at the V10 path."""
        if not (self.use_slice_context and self.use_context_adapter):
            return level
        rank = self.context_to_rank(context)
        adapter = torch.einsum("br,prh->bph", rank, self.context_adapter_basis)
        return level + torch.tanh(adapter)

    def forward(self, center, local, regional, local_position, regional_position,
                slice_descriptor, slice_mean, slice_std, local_key_padding_mask=None,
                regional_key_padding_mask=None, need_attention=False):
        batch = center.shape[0]
        if self.use_canonicalization:
            alpha = torch.sigmoid(self.canonical_logit)
            canonical_center = (center - slice_mean) / slice_std
            canonical_local = (local - slice_mean[:, None, :]) / slice_std[:, None, :]
            canonical_regional = (regional - slice_mean[:, None, :]) / slice_std[:, None, :]
            # A bounded residual update prevents an almost-constant gene from
            # dominating the encoder through a large inverse standard deviation.
            center = center + alpha * (canonical_center - center).clamp(-3.0, 3.0)
            local = local + alpha * (canonical_local - local).clamp(-3.0, 3.0)
            regional = regional + alpha * (canonical_regional - regional).clamp(-3.0, 3.0)
        context = self.context_encoder(slice_descriptor) if self.use_slice_context else torch.zeros(
            batch, self.center_norm.normalized_shape[0], device=center.device, dtype=center.dtype)
        z_center = self.center_norm(self.rna_encoder(center))                 # [B,H]
        local_tokens = self.encode_neighbours(local, local_position)           # [B,4,H]
        regional_tokens = self.encode_neighbours(regional, regional_position)  # [B,4,H]
        if self.protein_aware_attention_enabled:
            protein_id = torch.arange(self.n_proteins, device=center.device)
            protein_query = z_center[:, None, :] + self.protein_query_projection(
                self.protein_embedding(protein_id)
            )[None, :, :]                                                      # [B,P,H]
            z_local, local_attention = self.local_pool(
                protein_query, local_tokens, local_key_padding_mask, need_attention)
            z_regional, regional_attention = self.regional_pool(
                protein_query, regional_tokens, regional_key_padding_mask, need_attention)
        else:
            # A single centre-derived query is shared by all proteins.  This is
            # the direct generic-attention replacement for protein-aware queries.
            shared_query = z_center[:, None, :]
            z_local, local_attention = self.local_pool(
                shared_query, local_tokens, local_key_padding_mask, need_attention)
            z_regional, regional_attention = self.regional_pool(
                shared_query, regional_tokens, regional_key_padding_mask, need_attention)
            z_local = z_local.expand(-1, self.n_proteins, -1)
            z_regional = z_regional.expand(-1, self.n_proteins, -1)
        z_center_p = z_center[:, None, :].expand(batch, self.n_proteins, -1)
        levels = torch.stack((z_center_p, z_local, z_regional), dim=2)         # [B,P,3,H]
        router_input = torch.cat((z_center, z_local.mean(1), z_regional.mean(1)), dim=1)
        if self.use_slice_context:
            router_input = router_input + self.context_router(context)
        global_logits = self.global_router(router_input)[:, None, :]           # [B,1,3]
        if self.protein_specific_router_enabled:
            router_latent = self.router_context(router_input)                    # [B,R]
            protein_logits = torch.einsum("br,prs->bps", router_latent, self.protein_router)
            logits = global_logits + protein_logits + self.protein_router_bias[None]
        else:
            # Retain the global V9-style router, but remove all protein-specific
            # scale degrees of freedom.
            logits = global_logits.expand(-1, self.n_proteins, -1)
        weights = torch.softmax(logits, dim=-1)
        z_final = (levels * weights.unsqueeze(-1)).sum(dim=2)                  # [B,P,H]
        if self.use_cognate_anchor:
            u_c = center[:, self.cognate_indices]
            u_l = local[:, :, self.cognate_indices].mean(dim=1)
            u_r = regional[:, :, self.cognate_indices].mean(dim=1)
            u = torch.stack((u_c, u_l, u_r), dim=-1)

            z_f_ad = self.context_adapt(z_final, context)
            z_f_anc, d_f, gate = self.cognate_anchor(u, z_f_ad, self.has_cognate, scale_weights=weights)
            pred_final = self.protein_head(z_f_anc) + gate * d_f

            z_c_ad = self.context_adapt(z_center_p, context)
            z_c_anc, d_c, _ = self.cognate_anchor(u, z_c_ad, self.has_cognate)
            pred_center = self.protein_head(z_c_anc) + gate * d_c

            z_l_ad = self.context_adapt(z_local, context)
            z_l_anc, d_l, _ = self.cognate_anchor(u, z_l_ad, self.has_cognate)
            pred_local = self.protein_head(z_l_anc) + gate * d_l

            z_r_ad = self.context_adapt(z_regional, context)
            z_r_anc, d_r, _ = self.cognate_anchor(u, z_r_ad, self.has_cognate)
            pred_regional = self.protein_head(z_r_anc) + gate * d_r
        else:
            pred_final = self.protein_head(self.context_adapt(z_final, context))
            pred_center = self.protein_head(self.context_adapt(z_center_p, context))
            pred_local = self.protein_head(self.context_adapt(z_local, context))
            pred_regional = self.protein_head(self.context_adapt(z_regional, context))
            gate = None

        return {
            "final": pred_final,
            "center": pred_center,
            "local": pred_local,
            "regional": pred_regional,
            "protein_scale_attention": weights,
            "local_attention": local_attention,
            "regional_attention": regional_attention,
            "slice_context": context,
            "anchor_gate": gate,
        }


def protein_pcc_loss(prediction, target, eps=1e-6, min_variance=1e-6):
    """One minus mean protein-wise Pearson over a minibatch, with safeguards."""
    pred = prediction - prediction.mean(dim=0, keepdim=True)
    truth = target - target.mean(dim=0, keepdim=True)
    prediction_ss = pred.square().sum(0)
    target_ss = truth.square().sum(0)
    # Clamp *before* sqrt.  Clamping afterwards keeps the forward value finite
    # but leaves sqrt(0)'s singular derivative in the backward graph.  This
    # occurs on undersized final minibatches (e.g. two spots), where a protein
    # can have exactly zero batch variance.
    denominator = torch.sqrt((prediction_ss * target_ss).clamp_min(eps * eps))
    correlation = (pred * truth).sum(0) / denominator
    # Exclude a temporarily collapsed prediction head: its Pearson gradient is
    # ill-conditioned and otherwise can make the complete batch non-finite.
    valid = (torch.isfinite(correlation) & (target_ss > min_variance)
             & (prediction_ss > min_variance))
    if not valid.any():
        return prediction.new_zeros(()), 0
    return 1.0 - correlation[valid].mean(), int(valid.sum().item())


def protein_rmse_loss(prediction, target, eps=1e-8):
    """Mean batch RMSE over proteins; the same reduction as RMSE-Pro."""
    return torch.sqrt((prediction - target).square().mean(dim=0) + eps).mean()


def metric_rmse_pro(metric):
    return float(np.nanmean([row["rmse"] for row in metric["per_protein"]]))


def model_batch_inputs(x, local, regional, lp, rp, index, section_ids, slice_descriptors,
                       slice_gene_mean, slice_gene_std, device):
    """Prepare local V10 inputs plus RNA-only full-slice context for a batch."""
    from triscale_center_protein_v2_relpos import batch_inputs
    center, local_rna, regional_rna, local_pos, regional_pos = batch_inputs(
        x, local, regional, lp, rp, index, device)
    batch_sections = torch.as_tensor(section_ids[index], dtype=torch.long, device=device)
    return (center, local_rna, regional_rna, local_pos, regional_pos,
            slice_descriptors[batch_sections], slice_gene_mean[batch_sections],
            slice_gene_std[batch_sections])


def evaluate(model, x, y, local, regional, lp, rp, index, section_ids, slice_descriptors,
             slice_gene_mean, slice_gene_std, proteins, device, batch_size, details=False):
    from triscale_center_protein_v2_relpos import batch_inputs
    from triscale_center_protein_v9_dgat_compatible import dgat_metrics
    model.eval(); preds = []; scale = []; local_a = []; regional_a = []
    with torch.no_grad():
        for start in range(0, len(index), batch_size):
            idx = index[start:start + batch_size]
            out = model(*model_batch_inputs(
                x, local, regional, lp, rp, idx, section_ids, slice_descriptors,
                slice_gene_mean, slice_gene_std, device), need_attention=details)
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
    parser.add_argument("--pcc_min_variance", type=float, default=1e-6,
                        help="Minimum centered sum of squares for a protein to enter batch PCC loss.")
    parser.add_argument("--disable_protein_specific_router", action="store_true")
    parser.add_argument("--disable_protein_aware_attention", action="store_true")
    parser.add_argument("--disable_slice_context", action="store_true")
    parser.add_argument("--disable_canonicalization", action="store_true")
    parser.add_argument("--disable_context_adapter", action="store_true")
    parser.add_argument("--disable_cognate_anchor", action="store_true")
    parser.add_argument("--context_pca_components", type=int, default=32)
    parser.add_argument("--context_rank", type=int, default=8)
    parser.add_argument("--gene_dropout_rate", type=float, default=.5); parser.add_argument("--neighbor_mask_rate", type=float, default=.25)
    parser.add_argument("--local_neighbors", type=int, default=4); parser.add_argument("--regional_neighbors", type=int, default=8)
    parser.add_argument("--program_components", type=int, default=9); parser.add_argument("--save_epoch_checkpoints", action="store_true")
    args = parser.parse_args()
    if not 0 <= args.pcc_weight or args.pcc_min_variance <= 0 or args.context_pca_components < 2 or args.context_rank < 1 or not 0 <= args.gene_dropout_rate < 1 or not 0 <= args.neighbor_mask_rate < 1:
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
        if name in sources:
            adt = ad.read_h5ad(root / name / "adt_dgat_clr.h5ad")[:, proteins].copy()
            if not rna.obs_names.equals(adt.obs_names): raise ValueError(f"RNA/ADT barcode mismatch: {name}")
            block["y"] = dense(adt.X)
        blocks[name] = block
        audit[name] = {"spots": len(rna), "neighbour_sha256": digest(neighbours), "coordinate_scale": scale}
    target = blocks[args.target]; target_index = np.arange(len(target["x"]))
    parts = [[] for _ in range(8)]; section_parts = []; offset = 0; splits = {}
    for position, name in enumerate(sources):
        b = blocks[name]; order = np.arange(len(b["x"])); np.random.default_rng(args.seed + position).shuffle(order)
        n_val = max(1, int(round(len(order) * args.validation_ratio))); val, fit = np.sort(order[:n_val]), np.sort(order[n_val:])
        for dest, source in zip(parts[:6], [b["x"], b["y"], b["local"] + offset, b["regional"] + offset, b["lp"], b["rp"]]): dest.append(source)
        parts[6].append(fit + offset); parts[7].append(val + offset)
        section_parts.append(np.full(len(order), position, dtype=np.int64))
        splits[name] = {"fit": len(fit), "validation": len(val)}; offset += len(order)
    x, y, local, regional, lp, rp, fit, validation = map(np.concatenate, parts)
    section_ids = np.concatenate(section_parts)
    np.savez_compressed(out / "source_split_indices.npz", fit=fit, validation=validation)
    pca = PCA(n_components=args.program_components, svd_solver="full", random_state=args.seed).fit(y[fit])
    np.save(out / "pca_basis_source_fit.npy", pca.components_); np.save(out / "pca_mean_source_fit.npy", pca.mean_)

    # Fit RNA context axes on source fit spots only.  Each section descriptor
    # is a [mean, standard deviation] summary in this source-defined RNA space,
    # plus three RNA/geometry-only size descriptors.  No target protein is read.
    context_components = min(args.context_pca_components, len(fit) - 1, x.shape[1])
    context_pca = PCA(n_components=context_components, svd_solver="randomized", random_state=args.seed).fit(x[fit])
    np.save(out / "rna_slice_context_pca_basis_source_fit.npy", context_pca.components_)
    np.save(out / "rna_slice_context_pca_mean_source_fit.npy", context_pca.mean_)
    ordered_blocks = [blocks[name] for name in [*sources, args.target]]
    context_descriptors = []
    gene_mean = []
    gene_std = []
    for block in ordered_blocks:
        score = context_pca.transform(block["x"])
        coord = block["lp"][:, 0, :2]
        summary = np.array([
            np.log1p(len(block["x"])),
            np.log1p(np.linalg.norm(coord, axis=1).mean()),
            np.log1p(np.linalg.norm(coord, axis=1).std()),
        ], dtype=np.float32)
        context_descriptors.append(np.concatenate((score.mean(0), score.std(0), summary)).astype(np.float32))
        gene_mean.append(block["x"].mean(0).astype(np.float32))
        # The RNA matrix is already gene-scaled.  A conservative floor avoids
        # amplifying nearly constant genes during slice canonicalization.
        gene_std.append(np.clip(block["x"].std(0), 0.1, None).astype(np.float32))
    context_descriptors = np.stack(context_descriptors)
    gene_mean = np.stack(gene_mean)
    gene_std = np.stack(gene_std)
    # The descriptor encoder should see comparable dimensions but its fit uses
    # source sections only; applying the same transform to the target is valid.
    descriptor_center = context_descriptors[:-1].mean(0)
    # One-source transfer is a valid same-tissue setting.  Its descriptor
    # variance is exactly zero, so use a conservative source-only floor and a
    # bounded representation rather than allowing an unbounded target context.
    descriptor_scale = np.clip(context_descriptors[:-1].std(0), 0.1, None)
    context_descriptors = np.clip((context_descriptors - descriptor_center) / descriptor_scale, -5.0, 5.0)
    np.savez_compressed(out / "rna_slice_context_descriptors.npz", descriptors=context_descriptors,
                        source_center=descriptor_center, source_scale=descriptor_scale,
                        section_order=np.asarray([*sources, args.target]))
    np.savez_compressed(out / "rna_slice_canonicalization.npz", mean=gene_mean, std=gene_std,
                        section_order=np.asarray([*sources, args.target]))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    slice_descriptors_t = torch.from_numpy(context_descriptors).to(device)
    gene_mean_t = torch.from_numpy(gene_mean).to(device)
    gene_std_t = torch.from_numpy(gene_std).to(device)
    target_section_ids = np.full(len(target["x"]), len(sources), dtype=np.int64)
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
    cognate_indices = []
    has_cognate = []
    for p in proteins:
        g = cognate_map_dict.get(p)
        if g is not None and g in gene_to_idx:
            cognate_indices.append(gene_to_idx[g])
            has_cognate.append(1.0)
        else:
            cognate_indices.append(0)
            has_cognate.append(0.0)


    model = SliceConditionedNicheProtein(
        len(genes), len(proteins), args.hidden_dim, args.heads, args.ffn_dim,
        args.dropout, args.protein_embedding_dim, args.router_rank,
        context_input_dim=context_descriptors.shape[1], context_rank=args.context_rank,
        protein_specific_router=not args.disable_protein_specific_router,
        protein_aware_attention=not args.disable_protein_aware_attention,
        use_slice_context=not args.disable_slice_context,
        use_canonicalization=not args.disable_canonicalization,
        use_context_adapter=not args.disable_context_adapter,
        cognate_indices=cognate_indices,
        has_cognate=has_cognate,
        use_cognate_anchor=not args.disable_cognate_anchor,
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    # This diagnostic intentionally opens target ADT before training so that
    # every epoch can be scored on the target. It is an oracle/test-set
    # selection experiment and must never be used as a formal benchmark.
    target_adt = ad.read_h5ad(root / args.target / "adt_dgat_clr.h5ad")[:, proteins].copy()
    if not np.array_equal(np.asarray(target_adt.obs_names).astype(str), np.asarray(target["obs"]).astype(str)):
        raise ValueError(f"RNA/ADT barcode mismatch: {args.target}")
    target["y"] = dense(target_adt.X)
    loader = DataLoader(TensorDataset(torch.from_numpy(fit)), batch_size=args.batch_size, shuffle=True,
                        generator=torch.Generator().manual_seed(args.seed))
    best_val, best_epoch, best_state = -np.inf, 0, None; history = []
    oracle_best = {
        "pcc_pro_mean": {"value": -np.inf, "epoch": 0},
        "rmse_global": {"value": np.inf, "epoch": 0},
        "mae_global": {"value": np.inf, "epoch": 0},
    }
    oracle_best_metrics = {}
    print(f"SCN-Pro parameters={sum(v.numel() for v in model.parameters())}; sources={sources}; pcc_weight={args.pcc_weight}; context_pca={context_components}", flush=True)
    for epoch in range(1, args.epochs + 1):
        model.train(); total = []; final_losses = []; protein_rmse_losses = []; aux_losses = []; corr_losses = []; valid_counts = []
        for (idx_tensor,) in loader:
            idx = idx_tensor.numpy(); truth = torch.from_numpy(y[idx]).to(device)
            center, local_rna, regional_rna, local_pos, regional_pos, descriptor, slice_mean, slice_std = model_batch_inputs(
                x, local, regional, lp, rp, idx, section_ids, slice_descriptors_t, gene_mean_t, gene_std_t, device)
            output = model(feature_dropout(center, args.gene_dropout_rate), feature_dropout(local_rna, args.gene_dropout_rate),
                           feature_dropout(regional_rna, args.gene_dropout_rate), local_pos, regional_pos,
                           descriptor, slice_mean, slice_std,
                           sample_key_padding_mask(len(idx), local_rna.shape[1], args.neighbor_mask_rate, device),
                           sample_key_padding_mask(len(idx), regional_rna.shape[1], args.neighbor_mask_rate, device))
            final_mse = F.mse_loss(output["final"], truth)
            final_protein_rmse = protein_rmse_loss(output["final"], truth)
            auxiliary = sum(F.mse_loss(output[key], truth) for key in ("center", "local", "regional")) / 3
            corr, valid = protein_pcc_loss(output["final"], truth, min_variance=args.pcc_min_variance)
            loss = final_mse + args.aux_weight * auxiliary + args.pcc_weight * corr
            optimizer.zero_grad(set_to_none=True); loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step()
            total.append(loss.item()); final_losses.append(final_mse.item()); protein_rmse_losses.append(final_protein_rmse.item()); aux_losses.append(auxiliary.item()); corr_losses.append(corr.item()); valid_counts.append(valid)
        val, _, _ = evaluate(model, x, y, local, regional, lp, rp, validation, section_ids, slice_descriptors_t,
                              gene_mean_t, gene_std_t, proteins, device, args.batch_size)
        if args.save_epoch_checkpoints: torch.save(model.state_dict(), out / "epoch_checkpoints" / f"epoch_{epoch:02d}.pth")
        target_metric, target_pred, _ = evaluate(
            model, target["x"], target["y"], target["local"], target["regional"], target["lp"], target["rp"],
            target_index, target_section_ids, slice_descriptors_t, gene_mean_t, gene_std_t, proteins, device,
            args.batch_size)
        for metric_name, direction in (("pcc_pro_mean", "max"), ("rmse_global", "min"), ("mae_global", "min")):
            value = float(target_metric[metric_name])
            old = oracle_best[metric_name]["value"]
            improved = value > old + 1e-12 if direction == "max" else value < old - 1e-12
            if improved:
                oracle_best[metric_name] = {"value": value, "epoch": epoch}
                oracle_best_metrics[metric_name] = {k: (float(v) if isinstance(v, (np.floating, float)) else v)
                                                    for k, v in target_metric.items() if k != "per_protein"}
                np.save(out / f"target_oracle_best_{metric_name}_predictions.npy", target_pred)
        selection_score = val["pcc_pro_mean"]
        source_best = selection_score > best_val + 1e-6
        if source_best:
            best_val, best_epoch = selection_score, epoch
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            torch.save(best_state, out / "best_source_validation.pth")
        row = {"epoch": epoch, "train_loss": float(np.mean(total)), "train_final_mse": float(np.mean(final_losses)), "train_final_protein_rmse": float(np.mean(protein_rmse_losses)),
               "train_aux_mse": float(np.mean(aux_losses)), "train_pcc_loss": float(np.mean(corr_losses)), "valid_proteins_per_batch": float(np.mean(valid_counts)),
               "source_validation_pcc_pro_mean": val["pcc_pro_mean"], "source_validation_rmse_pro": metric_rmse_pro(val),
               "selection_score": selection_score, "is_selected_best": source_best,
               "target_epoch_pcc_pro_mean": target_metric["pcc_pro_mean"],
               "target_epoch_rmse_global": target_metric["rmse_global"],
               "target_epoch_mae_global": target_metric["mae_global"]}
        history.append(row); pd.DataFrame(history).to_csv(out / "metrics_history.csv", index=False)
        print(f"Epoch {epoch:02d}/{args.epochs} loss={row['train_loss']:.4f} source_val_PCC={selection_score:.4f}" + (" [SELECTED_BEST]" if source_best else ""), flush=True)
    model.load_state_dict(best_state); strict, strict_pred, strict_details = evaluate(
        model, target["x"], target["y"], target["local"], target["regional"], target["lp"], target["rp"],
        target_index, target_section_ids, slice_descriptors_t, gene_mean_t, gene_std_t, proteins, device,
        args.batch_size, details=True)
    strict_rmse_pro = metric_rmse_pro(strict)
    for prefix, metric, prediction, detail in [("source_selected", strict, strict_pred, strict_details)]:
        pd.DataFrame(metric.pop("per_protein")).to_csv(out / f"target_{prefix}_per_protein_metrics.csv", index=False); np.save(out / f"target_{prefix}_predictions.npy", prediction); np.savez_compressed(out / f"target_{prefix}_attention.npz", **detail)
    result = {"model": "SCN-Pro", "target": args.target, "sources": sources, "seed": args.seed, "epochs": args.epochs,
              "architecture": {"protein_specific_router": not args.disable_protein_specific_router, "protein_aware_neighbour_attention": not args.disable_protein_aware_attention, "protein_specific_head": f"vectorized {len(proteins)} independent nonlinear heads", "protein_embedding_dim": args.protein_embedding_dim, "router_rank": args.router_rank, "slice_context": not args.disable_slice_context, "niche_canonicalization": not args.disable_canonicalization, "context_adapter": not args.disable_context_adapter, "context_pca_components": context_components, "context_rank": args.context_rank},
              "loss": {"final_mse": 1.0, "protein_rmse": 0.0, "protein_pcc": args.pcc_weight, "auxiliary_mse": args.aux_weight},
              "selection": "source validation selects the formal checkpoint; target ADT is intentionally used for per-epoch oracle diagnostics", "checkpoint_selection": "source_validation", "preprocessing": "frozen DGAT registry; Muon CLR axis=0", "common_panel": {"genes": len(genes), "proteins": proteins}, "source_split": splits, "audit": audit,
              "slice_context_audit": {"source_fitted_rna_pca": True, "target_protein_used_for_checkpoint_selection": False, "target_protein_used_for_oracle_epoch_scoring": True, "target_full_slice_rna_used": True, "canonicalization_alpha": float(torch.sigmoid(model.canonical_logit).detach().cpu())}, "n_parameters": sum(v.numel() for v in model.parameters()), "best_epoch": best_epoch, "best_selection_pcc_pro_mean": best_val, "best_source_validation_pcc_pro_mean": best_val, "target_source_selected": strict, "target_source_selected_rmse_pro": strict_rmse_pro,
              "test_set_oracle": {"warning": "Target ADT was used after every epoch to select diagnostic checkpoints; this is not a valid formal benchmark.", "best_by_metric": oracle_best, "best_metrics": oracle_best_metrics},
              "wall_seconds": time.monotonic()}
    (out / "results.json").write_text(json.dumps(result, indent=2)); print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
