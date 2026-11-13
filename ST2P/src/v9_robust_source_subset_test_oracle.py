#!/usr/bin/env python3
"""Robust V9 source-subset diagnostic with RNA feature dropout and token masking.

The protected V9 architecture is instantiated unchanged. During training only,
RNA features are randomly dropped and selected neighbour tokens are excluded
from cross-attention. Every epoch is evaluated on held-out A1 for a clearly
labelled test-oracle diagnostic; A1 labels never enter optimisation.
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
from sklearn.neighbors import NearestNeighbors
from torch import nn
from torch.utils.data import DataLoader, TensorDataset, WeightedRandomSampler


def digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def corrected_knn(coords, k):
    """Query explicitly, then exclude self by identity (also handles ties)."""
    coords = np.asarray(coords)
    if not 0 < k < len(coords):
        raise ValueError(f"Require 0 < k < n_spots, got {k}, {len(coords)}")
    candidates = NearestNeighbors(n_neighbors=k + 1, n_jobs=1).fit(coords).kneighbors(
        coords, return_distance=False)
    index = np.stack([row[row != i][:k] for i, row in enumerate(candidates)])
    if index.shape != (len(coords), k) or np.any(index == np.arange(len(coords))[:, None]):
        raise AssertionError("Incorrect neighbour count or self retained")
    if not all(len(np.unique(row)) == k for row in index):
        raise AssertionError("Duplicate neighbour IDs")
    return index.astype(np.int64)


class IndependentProteinMLP(nn.Module):
    """Vectorized independent Linear(H,128)->LeakyReLU->Linear(128,1) heads.

    Each protein owns all its head weights. No BatchNorm, extra normalization,
    or dropout is added, so this factor changes the output parameterization.
    The return tuple matches ProgramResidualDecoder's interface.
    """
    def __init__(self, hidden, n_protein, width=128):
        super().__init__()
        self.weight1 = nn.Parameter(torch.empty(n_protein, hidden, width))
        self.bias1 = nn.Parameter(torch.empty(n_protein, width))
        self.weight2 = nn.Parameter(torch.empty(n_protein, width))
        self.bias2 = nn.Parameter(torch.empty(n_protein))
        # Match the per-Linear PyTorch default initialization.
        nn.init.uniform_(self.weight1, -1 / math.sqrt(hidden), 1 / math.sqrt(hidden))
        nn.init.uniform_(self.bias1, -1 / math.sqrt(hidden), 1 / math.sqrt(hidden))
        nn.init.uniform_(self.weight2, -1 / math.sqrt(width), 1 / math.sqrt(width))
        nn.init.uniform_(self.bias2, -1 / math.sqrt(width), 1 / math.sqrt(width))

    def forward(self, z):
        h = torch.nn.functional.leaky_relu(torch.einsum("bh,phw->bpw", z, self.weight1) + self.bias1)
        return (h * self.weight2).sum(-1) + self.bias2, None, None


def neighbour_pool_with_mask(pool, query, context, key_padding_mask, need_attention=False):
    """Run the existing V9 attention pool while excluding masked neighbours.

    ``key_padding_mask=True`` means a token is invisible to all attention
    heads. This reuses every protected V9 parameter and operation other than
    the optional masking argument supplied to ``MultiheadAttention``.
    """
    message, weight = pool.attention(
        pool.query_norm(query), pool.context_norm(context), pool.context_norm(context),
        key_padding_mask=key_padding_mask, need_weights=need_attention,
        average_attn_weights=False,
    )
    output = pool.norm_attention(message)
    output = pool.norm_ffn(output + pool.dropout(pool.ffn(output)))
    return output, weight


def robust_v9_forward(model, center, local, regional, local_position, regional_position,
                      local_key_padding_mask=None, regional_key_padding_mask=None,
                      need_attention=False):
    """Exact V9 forward path, with optional masked local/regional tokens."""
    z_center = model.center_norm(model.rna_encoder(center))
    query = z_center.unsqueeze(1)
    local_tokens = model.encode_neighbours(local, local_position)
    regional_tokens = model.encode_neighbours(regional, regional_position)
    z_local, local_attention = neighbour_pool_with_mask(
        model.local_pool, query, local_tokens, local_key_padding_mask, need_attention)
    z_regional, regional_attention = neighbour_pool_with_mask(
        model.regional_pool, query, regional_tokens, regional_key_padding_mask, need_attention)
    z_local, z_regional = z_local.squeeze(1), z_regional.squeeze(1)
    levels = torch.stack((z_center, z_local, z_regional), dim=1)
    level_weights = model.level_router(torch.cat((z_center, z_local, z_regional), dim=1))
    z_final = (levels * level_weights.unsqueeze(-1)).sum(dim=1)
    final, activity, residual = model.protein_program_decoder(z_final)
    center_pred, _, _ = model.protein_program_decoder(z_center)
    local_pred, _, _ = model.protein_program_decoder(z_local)
    regional_pred, _, _ = model.protein_program_decoder(z_regional)
    return {
        "final": final, "center": center_pred, "local": local_pred, "regional": regional_pred,
        "program_activity": activity, "residual": residual,
        "local_attention": local_attention, "regional_attention": regional_attention,
        "scale_attention": level_weights[:, None, None, :],
    }


def feature_dropout(x, rate):
    """Independent RNA feature perturbation, active during training only."""
    return F.dropout(x, p=rate, training=True) if rate > 0 else x


def sample_key_padding_mask(batch_size, tokens, rate, device):
    """Mask neighbour tokens but always retain at least one token per spot."""
    if rate <= 0:
        return None
    mask = torch.rand(batch_size, tokens, device=device) < rate
    all_masked = mask.all(dim=1)
    if all_masked.any():
        rows = torch.where(all_masked)[0]
        keep = torch.randint(tokens, (len(rows),), device=device)
        mask[rows, keep] = False
    return mask


def panel(root, names, filename):
    axes = []
    for name in names:
        a = ad.read_h5ad(root / name / filename, backed="r")
        axes.append(list(map(str, a.var_names)))
        a.file.close()
    common = set.intersection(*(set(v) for v in axes))
    return [v for v in axes[0] if v in common]


def dense(x):
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x, dtype=np.float32)


def prepare(root, names, target, genes, proteins, local_k, outer_k):
    from triscale_center_protein_v2_relpos import relative_position
    blocks, audits = {}, {}
    for name in names:
        r = ad.read_h5ad(root / name / "rna_dgat.h5ad")[:, genes].copy()
        coords = np.asarray(r.obsm["spatial"], dtype=np.float32)
        neighbours = corrected_knn(coords, outer_k)
        old = NearestNeighbors(n_neighbors=outer_k + 1).fit(coords).kneighbors(return_distance=False)[:, 1:]
        local, regional = neighbours[:, :local_k], neighbours[:, local_k:]
        distances = np.linalg.norm(coords[local[:, 0]] - coords, axis=1)
        scale = float(np.median(distances[distances > 0]))
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError(f"Invalid coordinate scale for {name}")
        block = dict(x=dense(r.X), local=local, regional=regional,
                     lp=relative_position(coords, local, scale),
                     rp=relative_position(coords, regional, scale), obs_names=np.asarray(r.obs_names))
        if name != target:
            a = ad.read_h5ad(root / name / "adt_dgat_clr.h5ad")[:, proteins].copy()
            if not r.obs_names.equals(a.obs_names):
                raise ValueError(f"RNA/ADT barcode order differs for {name}")
            block["y"] = dense(a.X)
        blocks[name] = block
        audits[name] = dict(spots=len(r), coordinate_scale=scale, neighbours_sha256=digest(neighbours),
                            self_count=int((neighbours == np.arange(len(r))[:, None]).sum()),
                            local_count=local_k, regional_count=outer_k - local_k,
                            fraction_spots_neighbour_set_changed=float(np.mean([
                                set(a) != set(b) for a, b in zip(neighbours, old)])))
        print(f"Prepared {name}: {len(r)} spots, corrected local={local_k}/regional={outer_k-local_k}", flush=True)
    return blocks, audits


def evaluate(model, x, y, local, regional, lp, rp, index, proteins, device, batch_size):
    from triscale_center_protein_v2_relpos import batch_inputs
    from triscale_center_protein_v9_dgat_compatible import dgat_metrics
    model.eval()
    prediction, scales = [], []
    with torch.no_grad():
        for start in range(0, len(index), batch_size):
            idx = index[start:start + batch_size]
            output = robust_v9_forward(model, *batch_inputs(x, local, regional, lp, rp, idx, device))
            prediction.append(output["final"].cpu().numpy())
            scales.append(output["scale_attention"].cpu().numpy())
    pred = np.concatenate(prediction)
    return dgat_metrics(pred, y[index], proteins), pred, np.concatenate(scales)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--prepared_root", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--experiment_dir", required=True)
    p.add_argument("--head", choices=["program", "per_protein_mlp"], required=True)
    p.add_argument("--aux_weight", type=float, required=True)
    p.add_argument("--gene_dropout_rate", type=float, default=.5,
                   help="Training-only RNA feature dropout; .5 matches NicheTrans's noise dropout.")
    p.add_argument("--neighbor_mask_rate", type=float, default=.25,
                   help="Training-only probability of masking each neighbour token.")
    p.add_argument("--target", default="A1_LN")
    p.add_argument("--source_sections", nargs="+", default=None,
                   help="Explicit source sections. Omit to use every non-target section.")
    p.add_argument("--d1_sampling_probability", type=float, default=None,
                   help=("Optional training-only D1_LN mixture weight. For example, .5 samples "
                         "D1_LN with probability 50%% and all other source spots collectively with "
                         "probability 50%%. Source validation stays unweighted and section-stratified."))
    p.add_argument("--save_epoch_checkpoints", action="store_true",
                   help="Save model weights after every epoch for retrospective diagnostics.")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--epochs", type=int, default=20)
    p.add_argument("--batch_size", type=int, default=64)
    p.add_argument("--validation_ratio", type=float, default=.1)
    p.add_argument("--hidden_dim", type=int, default=256)
    p.add_argument("--heads", type=int, default=4)
    p.add_argument("--ffn_dim", type=int, default=256)
    p.add_argument("--dropout", type=float, default=.3)
    p.add_argument("--lr", type=float, default=5e-4)
    p.add_argument("--weight_decay", type=float, default=1e-4)
    p.add_argument("--program_components", type=int, default=9)
    p.add_argument("--local_neighbors", type=int, default=6)
    p.add_argument("--regional_neighbors", type=int, default=24)
    args = p.parse_args()
    if not 0 <= args.gene_dropout_rate < 1 or not 0 <= args.neighbor_mask_rate < 1:
        raise ValueError("Dropout and neighbour-mask rates must be in [0, 1).")
    sys.path.insert(0, str(Path(args.experiment_dir).resolve()))
    from triscale_center_protein_v9_dgat_compatible import TriScaleCenterProteinV9
    from triscale_center_protein_v2_relpos import batch_inputs

    random.seed(args.seed); np.random.seed(args.seed); torch.manual_seed(args.seed)
    if torch.cuda.is_available(): torch.cuda.manual_seed_all(args.seed)
    torch.backends.cudnn.deterministic = True; torch.backends.cudnn.benchmark = False
    start_time = time.monotonic()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=False)
    epoch_checkpoint_dir = out / "epoch_checkpoints"
    if args.save_epoch_checkpoints:
        epoch_checkpoint_dir.mkdir()
    (out / "config.json").write_text(json.dumps(vars(args), indent=2))
    shutil.copy2(__file__, out / "script_snapshot.py")
    root = Path(args.prepared_root)
    manifest = json.loads((root / "prepared_manifest.json").read_text())
    if manifest.get("adt_clr_audit", {}).get("axis") != 0:
        raise ValueError("Require the frozen exact Muon axis=0 registry")
    names = [d["dataset"] for d in manifest["datasets"]]
    available_sources = [n for n in names if n != args.target]
    if args.source_sections is None:
        sources = available_sources
    else:
        sources = list(args.source_sections)
        unknown = sorted(set(sources) - set(available_sources))
        if unknown:
            raise ValueError(f"Invalid/non-source requested: {unknown}")
        if len(set(sources)) != len(sources):
            raise ValueError("Source sections must be unique")
    if args.target not in names: raise ValueError("Target absent")
    if args.d1_sampling_probability is not None:
        if not 0 < args.d1_sampling_probability < 1:
            raise ValueError("--d1_sampling_probability must lie strictly between 0 and 1")
        if "D1_LN" not in sources:
            raise ValueError("D1-weighted sampling requires D1_LN among --source_sections")
        if len(sources) < 2:
            raise ValueError("D1-weighted sampling requires at least one non-D1 source")
    genes = panel(root, names, "rna_dgat.h5ad")
    proteins = panel(root, names, "adt_dgat_clr.h5ad")
    blocks, audits = prepare(root, names, args.target, genes, proteins,
                             args.local_neighbors, args.regional_neighbors)
    # Deliberate test-oracle diagnostic: target labels are used for evaluation
    # after every epoch, but never enter gradients, optimizer updates, or loss.
    tb = blocks[args.target]
    target_adt = ad.read_h5ad(root / args.target / "adt_dgat_clr.h5ad")[:, proteins].copy()
    assert np.array_equal(tb["obs_names"], np.asarray(target_adt.obs_names))
    target_y = dense(target_adt.X)
    target_index = np.arange(len(tb["x"]))
    x_parts, y_parts, l_parts, r_parts, lp_parts, rp_parts, fits, vals = ([] for _ in range(8))
    fit_section_parts = []
    offset, splits = 0, {}
    for position, name in enumerate(sources):
        b = blocks[name]; order = np.arange(len(b["x"]))
        np.random.default_rng(args.seed + position).shuffle(order)
        n_val = max(1, int(round(len(order) * args.validation_ratio)))
        fit, val = np.sort(order[n_val:]), np.sort(order[:n_val])
        x_parts.append(b["x"]); y_parts.append(b["y"])
        l_parts.append(b["local"] + offset); r_parts.append(b["regional"] + offset)
        lp_parts.append(b["lp"]); rp_parts.append(b["rp"])
        fits.append(fit + offset); vals.append(val + offset)
        fit_section_parts.append(np.full(len(fit), name, dtype=object))
        splits[name] = dict(fit=len(fit), validation=len(val))
        offset += len(order)
    x, y, local, regional, lp, rp, fit, val_idx = map(np.concatenate,
        (x_parts, y_parts, l_parts, r_parts, lp_parts, rp_parts, fits, vals))
    fit_sections = np.concatenate(fit_section_parts)
    assert not np.intersect1d(fit, val_idx).size
    legacy = root.parent.parent / "results/dgat_muon017_axis0_a1_lodo_20261020/v9_seed42_official20epoch_budget/source_split_indices.npz"
    if args.source_sections is None and legacy.exists():
        with np.load(legacy) as old:
            assert np.array_equal(fit, old["fit"]) and np.array_equal(val_idx, old["validation"])
    np.savez_compressed(out / "source_split_indices.npz", fit=fit, validation=val_idx)
    pca = PCA(n_components=args.program_components, svd_solver="full", random_state=args.seed).fit(y[fit])
    np.save(out / "pca_basis_source_fit.npy", pca.components_)
    np.save(out / "pca_mean_source_fit.npy", pca.mean_)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = TriScaleCenterProteinV9(len(genes), pca.components_.astype(np.float32), pca.mean_.astype(np.float32),
                                  args.hidden_dim, args.heads, args.ffn_dim, args.dropout)
    encoder_hash = digest(torch.cat([q.detach().reshape(-1) for n, q in model.named_parameters()
                                     if not n.startswith("protein_program_decoder.")]).numpy())
    if args.head == "per_protein_mlp":
        model.protein_program_decoder = IndependentProteinMLP(args.hidden_dim, len(proteins))
    model.to(device)

    # With no perturbation/mask, this explicit forward must equal protected V9.
    model.eval()
    with torch.no_grad():
        smoke_idx = fit[:min(2, len(fit))]
        smoke_inputs = batch_inputs(x, local, regional, lp, rp, smoke_idx, device)
        reference = model(*smoke_inputs)
        masked_free = robust_v9_forward(model, *smoke_inputs)
        for key in ("final", "center", "local", "regional", "scale_attention"):
            if not torch.allclose(reference[key], masked_free[key], rtol=1e-6, atol=1e-6):
                raise AssertionError(f"Masked-forward differs from V9 without masks: {key}")
    model.train()
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    # Optional section-balanced sampling affects only optimization batches.  It
    # intentionally does not alter section-stratified source validation or any
    # A1 diagnostic evaluation.  Each epoch still performs len(fit) draws.
    sampling_audit = {"mode": "uniform_spot", "draws_per_epoch": int(len(fit))}
    if args.d1_sampling_probability is None:
        loader = DataLoader(TensorDataset(torch.from_numpy(fit)), batch_size=args.batch_size, shuffle=True,
                            generator=torch.Generator().manual_seed(args.seed))
    else:
        is_d1 = fit_sections == "D1_LN"
        n_d1, n_other = int(is_d1.sum()), int((~is_d1).sum())
        if n_d1 == 0 or n_other == 0:
            raise AssertionError("Missing D1 or non-D1 fit spots for weighted sampling")
        # Per-spot probabilities sum to the requested section mass, preventing
        # larger non-LN sections from dominating the source mixture.
        draw_weights = np.where(
            is_d1,
            args.d1_sampling_probability / n_d1,
            (1.0 - args.d1_sampling_probability) / n_other,
        ).astype(np.float64)
        sampler = WeightedRandomSampler(
            torch.from_numpy(draw_weights), num_samples=len(fit), replacement=True,
            generator=torch.Generator().manual_seed(args.seed),
        )
        loader = DataLoader(TensorDataset(torch.from_numpy(fit)), batch_size=args.batch_size, sampler=sampler)
        sampling_audit = {
            "mode": "D1_LN_weighted_with_replacement",
            "D1_LN_probability": args.d1_sampling_probability,
            "all_other_sources_probability": 1.0 - args.d1_sampling_probability,
            "D1_LN_fit_spots": n_d1,
            "other_fit_spots": n_other,
            "draws_per_epoch": int(len(fit)),
        }
    history = []
    best_epoch, best_pcc = 0, -np.inf
    best_test_epoch, best_test_pcc = 0, -np.inf
    print(f"START robust8 head={args.head} aux={args.aux_weight} gene_dropout={args.gene_dropout_rate} "
          f"neighbor_mask={args.neighbor_mask_rate} panel={len(genes)}/{len(proteins)} "
          f"parameters={sum(q.numel() for q in model.parameters())} encoder_init={encoder_hash}", flush=True)
    for epoch in range(1, args.epochs + 1):
        epoch_time = time.monotonic(); model.train()
        losses, final_losses, aux_losses = [], [], []
        for (ix,) in loader:
            idx = ix.numpy(); truth = torch.from_numpy(y[idx]).to(device)
            center, local_rna, regional_rna, local_pos, regional_pos = batch_inputs(
                x, local, regional, lp, rp, idx, device)
            output = robust_v9_forward(
                model,
                feature_dropout(center, args.gene_dropout_rate),
                feature_dropout(local_rna, args.gene_dropout_rate),
                feature_dropout(regional_rna, args.gene_dropout_rate),
                local_pos, regional_pos,
                sample_key_padding_mask(len(idx), local_rna.shape[1], args.neighbor_mask_rate, device),
                sample_key_padding_mask(len(idx), regional_rna.shape[1], args.neighbor_mask_rate, device),
            )
            final_loss = nn.functional.mse_loss(output["final"], truth)
            aux = sum(nn.functional.mse_loss(output[k], truth) for k in ("center", "local", "regional")) / 3
            loss = final_loss + args.aux_weight * aux
            optimizer.zero_grad(set_to_none=True); loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step()
            losses.append(loss.item()); final_losses.append(final_loss.item()); aux_losses.append(aux.item())
        val, _, _ = evaluate(model, x, y, local, regional, lp, rp, val_idx, proteins, device, args.batch_size)
        test, _, _ = evaluate(model, tb["x"], target_y, tb["local"], tb["regional"],
                              tb["lp"], tb["rp"], target_index, proteins, device, args.batch_size)
        if args.save_epoch_checkpoints:
            torch.save({k: v.detach().cpu().clone() for k, v in model.state_dict().items()},
                       epoch_checkpoint_dir / f"epoch_{epoch:02d}.pth")
        improved = val["pcc_pro_mean"] > best_pcc + 1e-6
        if improved:
            best_pcc, best_epoch = val["pcc_pro_mean"], epoch
            torch.save({k: v.detach().cpu().clone() for k, v in model.state_dict().items()},
                       out / "best_source_validation.pth")
            (out / "best_source_validation_metrics.json").write_text(json.dumps(val, indent=2))
        test_improved = test["pcc_pro_mean"] > best_test_pcc + 1e-6
        if test_improved:
            best_test_pcc, best_test_epoch = test["pcc_pro_mean"], epoch
            torch.save({k: v.detach().cpu().clone() for k, v in model.state_dict().items()},
                       out / "best_test_oracle.pth")
            (out / "best_test_oracle_metrics.json").write_text(json.dumps(test, indent=2))
        row = dict(epoch=epoch, train_loss=float(np.mean(losses)), train_final_mse=float(np.mean(final_losses)),
                   train_aux_mse=float(np.mean(aux_losses)), is_source_best=improved,
                   is_test_oracle_best=test_improved, seconds=time.monotonic()-epoch_time)
        row.update({"source_validation_"+k: v for k,v in val.items() if k != "per_protein"})
        row.update({"A1_test_"+k: v for k,v in test.items() if k != "per_protein"})
        history.append(row); pd.DataFrame(history).to_csv(out / "metrics_history.csv", index=False)
        print(f"Epoch {epoch:02d}/{args.epochs} loss={row['train_loss']:.5f} "
              f"source_val_PCCPro={val['pcc_pro_mean']:.4f} "
              f"A1_test_PCCPro={test['pcc_pro_mean']:.4f} "
              f"A1_RMSE={test['rmse_global']:.4f} seconds={row['seconds']:.1f}" +
              (" [SOURCE_BEST]" if improved else "") +
              (" [TEST_ORACLE_BEST]" if test_improved else ""), flush=True)

    # Save both protocols separately: source-selected remains strict; test-best
    # is explicitly oracle-only and is never mixed with the strict result.
    model.load_state_dict(torch.load(out / "best_source_validation.pth", map_location=device, weights_only=True))
    strict_metrics, strict_prediction, strict_scales = evaluate(
        model, tb["x"], target_y, tb["local"], tb["regional"], tb["lp"], tb["rp"],
        target_index, proteins, device, args.batch_size)
    pd.DataFrame(strict_metrics.pop("per_protein")).to_csv(
        out / "a1_source_selected_per_protein_metrics.csv", index=False)
    np.save(out / "a1_source_selected_predictions.npy", strict_prediction)
    np.save(out / "a1_source_selected_scale_attention.npy", strict_scales)

    model.load_state_dict(torch.load(out / "best_test_oracle.pth", map_location=device, weights_only=True))
    oracle_metrics, oracle_prediction, oracle_scales = evaluate(
        model, tb["x"], target_y, tb["local"], tb["regional"], tb["lp"], tb["rp"],
        target_index, proteins, device, args.batch_size)
    pd.DataFrame(oracle_metrics.pop("per_protein")).to_csv(
        out / "a1_test_oracle_per_protein_metrics.csv", index=False)
    np.save(out / "a1_test_oracle_predictions.npy", oracle_prediction)
    np.save(out / "a1_test_oracle_scale_attention.npy", oracle_scales)
    dependency_names = ["triscale_center_protein_v9_dgat_compatible.py", "triscale_center_protein_v5_program_decoder.py",
                        "triscale_center_protein_v7_anchored_cross_attention.py", "triscale_center_protein_v2_relpos.py"]
    dependencies = {n: hashlib.sha256((Path(args.experiment_dir)/n).read_bytes()).hexdigest() for n in dependency_names}
    result = dict(model="V9 corrected-neighbour ablation", head=args.head, aux_weight=args.aux_weight,
                  target=args.target, sources=sources, seed=args.seed, epochs=args.epochs,
                  epoch_checkpoints_saved=bool(args.save_epoch_checkpoints),
                  regularization=dict(gene_dropout_rate=args.gene_dropout_rate,
                                      neighbor_mask_rate=args.neighbor_mask_rate,
                                      active="training only; at least one neighbour retained per scale"),
                  training_sampling=sampling_audit,
                  selection=dict(
                      strict="source-only validation mean PCC-Pro",
                      diagnostic="A1 test PCC-Pro inspected every epoch; test-oracle only"),
                  preprocessing="Frozen DGAT registry, exact Muon 0.1.7 CLR axis=0; no extra normalization",
                  common_panel=dict(genes=len(genes), proteins=proteins), source_split=splits,
                  section_audit=audits, encoder_init_sha256=encoder_hash, dependency_sha256=dependencies,
                  fit_indices_sha256=digest(fit), validation_indices_sha256=digest(val_idx),
                  n_parameters=sum(q.numel() for q in model.parameters()), best_epoch=best_epoch,
                  best_source_validation_pcc_pro_mean=best_pcc,
                  source_selected_test_A1=strict_metrics,
                  best_test_oracle_epoch=best_test_epoch,
                  best_test_oracle_pcc_pro_mean=best_test_pcc,
                  test_oracle_A1=oracle_metrics,
                  target_adt_usage="evaluation-only every epoch; never used in training loss",
                  wall_seconds=time.monotonic()-start_time)
    (out / "results.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
