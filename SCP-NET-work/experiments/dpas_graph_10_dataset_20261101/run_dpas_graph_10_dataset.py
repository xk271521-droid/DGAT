"""Run DPAS-Graph on the fixed ten-dataset benchmark.

The model/data modules are imported from the pinned official source. This file
only adapts the project registry and enforces source-only checkpoint selection.
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import random
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from torch_geometric.loader import DataLoader


def find_project_root() -> Path:
    here = Path(__file__).resolve()
    candidates = [here.parent, *here.parents]
    for candidate in candidates:
        if (candidate / "ST2P" / "data" / "processed").exists() and (candidate / "upstream").exists():
            return candidate
    # Local checkout fallback: E:\DGAT\...
    return here.parents[3]


PROJECT_ROOT = find_project_root()
UPSTREAM_ROOT = PROJECT_ROOT / "upstream" / "DPAS-Graph"
DATA_ROOT = PROJECT_ROOT / "ST2P" / "data" / "processed" / "DGAT_final10_common9384_muon017_axis0_20261101"
PANEL_PATH = PROJECT_ROOT / "SCP-NET-work" / "docs" / "ten_dataset_common_panel_20261101.json"

PACKAGE_ROOT = UPSTREAM_ROOT / "dpas-graph"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

# The pinned repository directory is named ``dpas-graph`` while its source
# imports the package as ``dpas``. Register that directory under the author's
# import name without copying or editing any upstream file.
import importlib.util

if "dpas" not in sys.modules:
    _dpas_spec = importlib.util.spec_from_file_location(
        "dpas",
        PACKAGE_ROOT / "__init__.py",
        submodule_search_locations=[str(PACKAGE_ROOT)],
    )
    _dpas_module = importlib.util.module_from_spec(_dpas_spec)
    sys.modules["dpas"] = _dpas_module
    _dpas_spec.loader.exec_module(_dpas_module)

from dpas.data.graph_dataset import MultiGraphDataset_mRNA_with_protein_target
from dpas.eval.metrics import evaluate_prediction_arrays, save_prediction_artifacts
from dpas.models.dpas_graph import Decoder_Protein_MLP, Decoder_mRNA_MaskedSmall, DualAdaptiveEncoder


LAMBDA_PRED = 8.0
LAMBDA_RNA = 1.0
HIDDEN_DIM = 1024
DROPOUT_RATE = 0.3
BATCH_SIZE = 1


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def as_dense_float32(x) -> np.ndarray:
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x, dtype=np.float32)


def load_panel() -> tuple[list[str], list[str]]:
    panel = json.loads(PANEL_PATH.read_text(encoding="utf-8"))
    genes = [str(x) for x in panel["rna_common_genes"]]
    proteins = [str(x) for x in panel["adt_common_proteins"]]
    if len(genes) != 9384 or len(proteins) != 31:
        raise RuntimeError(f"Unexpected fixed panel: genes={len(genes)}, proteins={len(proteins)}")
    return genes, proteins


def load_dataset(name: str, genes: list[str], proteins: list[str], include_adt: bool = True):
    base = DATA_ROOT / name
    rna = sc.read_h5ad(base / "rna_dgat.h5ad")
    adt = sc.read_h5ad(base / "adt_dgat_clr.h5ad") if include_adt else None
    if adt is None:
        rna = rna[:, genes].copy()
        rna.uns["name"] = name
        rna.uns["_dpas_input_domain"] = "proc"
        if "spatial" not in rna.obsm:
            raise RuntimeError(f"{name}: RNA is missing obsm['spatial']")
        rna.X = as_dense_float32(rna.X)
        return rna, None
    common_spots = rna.obs_names.intersection(adt.obs_names)
    if len(common_spots) == 0:
        raise RuntimeError(f"{name}: RNA/ADT have no common spots")
    rna = rna[common_spots, genes].copy()
    adt = adt[common_spots, proteins].copy()
    rna.uns["name"] = name
    adt.uns["name"] = name
    rna.uns["_dpas_input_domain"] = "proc"
    adt.uns["_dpas_input_domain"] = "proc"
    if "spatial" not in rna.obsm:
        raise RuntimeError(f"{name}: RNA is missing obsm['spatial']")
    if "spatial" not in adt.obsm:
        adt.obsm["spatial"] = rna.obsm["spatial"].copy()
    rna.X = as_dense_float32(rna.X)
    adt.X = as_dense_float32(adt.X)
    return rna, adt


def split_source(rna, adt, seed: int, fraction: float = 0.10):
    n = rna.n_obs
    if n < 20:
        raise RuntimeError(f"{rna.uns['name']}: too few spots for source validation ({n})")
    rng = np.random.default_rng(seed)
    indices = np.arange(n)
    rng.shuffle(indices)
    n_val = max(1, int(round(n * fraction)))
    val_idx = np.sort(indices[:n_val])
    fit_idx = np.sort(indices[n_val:])
    return rna[fit_idx].copy(), adt[fit_idx].copy(), rna[val_idx].copy(), adt[val_idx].copy()


def make_loader(rnas, adts, path: Path, shuffle: bool):
    dataset = MultiGraphDataset_mRNA_with_protein_target(rnas, adts, save_dir=str(path))
    return DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=shuffle)


def model_bundle(n_genes: int, n_proteins: int, device):
    enc = DualAdaptiveEncoder(in_channels=n_genes, hidden_dim=HIDDEN_DIM, dropout=DROPOUT_RATE).to(device)
    dec_rna = Decoder_mRNA_MaskedSmall(hidden_dim=HIDDEN_DIM, output_dim=n_genes).to(device)
    dec_prot = Decoder_Protein_MLP(hidden_dim=HIDDEN_DIM, proteins_or_dim=n_proteins).to(device)
    return enc, dec_rna, dec_prot


def forward_batch(data, enc, dec_rna, dec_prot, device):
    m_x = data["mRNA"].x.to(device)
    m_ei = data[("mRNA", "mRNA_knn", "mRNA")].edge_index.to(device)
    m_ea = data[("mRNA", "mRNA_knn", "mRNA")].edge_attr.to(device)
    target = data.protein_target.to(device)
    z = enc(m_x, m_ei, m_ea)
    return z, dec_rna(z), dec_prot(z), target, m_x


def predict(loader, enc, dec_prot, device):
    enc.eval()
    dec_prot.eval()
    preds, truths = [], []
    with torch.no_grad():
        for data in loader:
            m_x = data["mRNA"].x.to(device)
            m_ei = data[("mRNA", "mRNA_knn", "mRNA")].edge_index.to(device)
            m_ea = data[("mRNA", "mRNA_knn", "mRNA")].edge_attr.to(device)
            z = enc(m_x, m_ei, m_ea)
            preds.append(dec_prot(z).detach().cpu().numpy())
            truths.append(data.protein_target.detach().cpu().numpy())
    return np.vstack(preds).astype(np.float32), np.vstack(truths).astype(np.float32)


def pcc_pro_mean(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    result = evaluate_prediction_arrays(y_true, y_pred, include_spearman=False)
    return float(result["summary"]["pcc_pro_mean"])


def save_state(path: Path, enc, dec_rna, dec_prot):
    path.mkdir(parents=True, exist_ok=True)
    torch.save(enc.state_dict(), path / "best_encoder_mRNA.pth")
    torch.save(dec_rna.state_dict(), path / "best_decoder_mRNA.pth")
    torch.save(dec_prot.state_dict(), path / "best_decoder_protein.pth")


def train_direction(name: str, sources: list[str], targets: list[str], args, genes, proteins):
    out = Path(args.out_root) / name
    out.mkdir(parents=True, exist_ok=True)
    log_path = out / "run.log"
    def log(msg):
        print(msg, flush=True)
        with log_path.open("a", encoding="utf-8") as f:
            f.write(str(msg) + "\n")

    set_seed(args.seed)
    loaded = {n: load_dataset(n, genes, proteins, include_adt=True) for n in sources}
    target_rna_only = {n: load_dataset(n, genes, proteins, include_adt=False)[0] for n in targets}
    fit_rnas, fit_adts, val_rnas, val_adts = [], [], [], []
    for i, src in enumerate(sources):
        fr, fa, vr, va = split_source(*loaded[src], seed=args.seed + i)
        fit_rnas.append(fr); fit_adts.append(fa); val_rnas.append(vr); val_adts.append(va)

    graphs = out / "graphs"
    fit_loader = make_loader(fit_rnas, fit_adts, graphs / "fit", shuffle=True)
    val_loader = make_loader(val_rnas, val_adts, graphs / "val", shuffle=False)
    target_rnas = [target_rna_only[n] for n in targets]

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    enc, dec_rna, dec_prot = model_bundle(len(genes), len(proteins), device)
    opt_enc = torch.optim.Adam(enc.parameters(), lr=5e-4 * args.lr_scale, weight_decay=2e-5)
    opt_rna = torch.optim.Adam(dec_rna.parameters(), lr=1e-4 * args.lr_scale, weight_decay=1e-5)
    opt_prot = torch.optim.Adam(dec_prot.parameters(), lr=1e-4 * args.lr_scale, weight_decay=1e-5)
    huber = torch.nn.SmoothL1Loss(beta=1.0)
    best_score = -float("inf")
    best_epoch = 0
    best_state = None
    history = []
    log(f"[{name}] device={device} sources={sources} targets={targets} genes={len(genes)} proteins={len(proteins)}")

    for epoch in range(1, args.epochs + 1):
        enc.train(); dec_rna.train(); dec_prot.train()
        total_loss = 0.0; count = 0
        for data in fit_loader:
            _, rna_pred, prot_pred, truth, m_x = forward_batch(data, enc, dec_rna, dec_prot, device)
            loss = LAMBDA_PRED * huber(prot_pred, truth) + LAMBDA_RNA * torch.sqrt(torch.mean((rna_pred - m_x) ** 2) + 1e-8)
            opt_enc.zero_grad(set_to_none=True); opt_rna.zero_grad(set_to_none=True); opt_prot.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(enc.parameters(), 1.0)
            torch.nn.utils.clip_grad_norm_(dec_rna.parameters(), 1.0)
            torch.nn.utils.clip_grad_norm_(dec_prot.parameters(), 1.0)
            opt_enc.step(); opt_rna.step(); opt_prot.step()
            total_loss += float(loss.detach().cpu()); count += 1
        val_pred, val_true = predict(val_loader, enc, dec_prot, device)
        score = pcc_pro_mean(val_true, val_pred)
        row = {"epoch": epoch, "train_loss": total_loss / max(1, count), "source_validation_pcc_pro_mean": score}
        history.append(row)
        log(f"[{name}] epoch={epoch}/{args.epochs} loss={row['train_loss']:.5f} val_pcc={score:.4f}")
        if score > best_score:
            best_score = score; best_epoch = epoch
            best_state = {"enc": copy.deepcopy(enc.state_dict()), "rna": copy.deepcopy(dec_rna.state_dict()), "prot": copy.deepcopy(dec_prot.state_dict())}
            save_state(out, enc, dec_rna, dec_prot)
        if epoch % 10 == 0:
            for opt in (opt_enc, opt_rna, opt_prot):
                for group in opt.param_groups: group["lr"] *= 0.8

    enc.load_state_dict(best_state["enc"]); dec_rna.load_state_dict(best_state["rna"]); dec_prot.load_state_dict(best_state["prot"])
    # Target ADT is opened only after source-only checkpoint selection.
    target_adts = [load_dataset(n, genes, proteins, include_adt=True)[1] for n in targets]
    target_loaders = [make_loader([r], [a], graphs / f"target_{n}", shuffle=False) for n, r, a in zip(targets, target_rnas, target_adts)]
    target_metrics = {}
    for target_name, loader in zip(targets, target_loaders):
        pred, true = predict(loader, enc, dec_prot, device)
        spot_names = list(target_rna_only[target_name].obs_names)
        eval_dict = evaluate_prediction_arrays(true, pred, protein_names=proteins, spot_names=spot_names)
        save_prediction_artifacts(out / target_name, eval_dict, target_name, protein_names=proteins, spot_names=spot_names, summary_extra={"direction": name, "selection_metric": "source_validation_pcc_pro_mean", "best_epoch": best_epoch})
        target_metrics[target_name] = eval_dict["summary"]

    (out / "metrics_history.json").write_text(json.dumps(history, indent=2), encoding="utf-8")
    summary = {"name": name, "sources": sources, "targets": targets, "genes": len(genes), "proteins": len(proteins), "best_epoch": best_epoch, "best_source_validation_pcc_pro_mean": best_score, "target_metrics": target_metrics, "official_parameters": {"lambda_pred": LAMBDA_PRED, "lambda_rna": LAMBDA_RNA, "hidden_dim": HIDDEN_DIM, "dropout_rate": DROPOUT_RATE, "epochs": args.epochs, "lr_scale": args.lr_scale, "seed": args.seed}, "device": str(device), "selection_metric": "source_validation_pcc_pro_mean"}
    (out / "results.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    log(f"[{name}] COMPLETE best_epoch={best_epoch} targets={targets}")


DIRECTIONS = [
    ("same_LN_A1_to_D1", ["A1_LN"], ["D1_LN"]), ("same_LN_D1_to_A1", ["D1_LN"], ["A1_LN"]),
    ("same_TNSL_A1_to_D1", ["A1_TNSL"], ["D1_TNSL"]), ("same_TNSL_D1_to_A1", ["D1_TNSL"], ["A1_TNSL"]),
    ("same_HumanTonsil_10x_to_AddOns", ["HumanTonsil_10x"], ["HumanTonsil_AddOns_S1"]), ("same_HumanTonsil_AddOns_to_10x", ["HumanTonsil_AddOns_S1"], ["HumanTonsil_10x"]),
    ("same_PBC_6837_to_6835", ["PBC_PR_6837"], ["PBC_PR_6835_5A"]), ("same_PBC_6835_to_6837", ["PBC_PR_6835_5A"], ["PBC_PR_6837"]),
    ("cross_holdout_LN", ["A1_TNSL", "D1_TNSL", "HumanTonsil_10x", "HumanTonsil_AddOns_S1", "Breast_Cancer", "PBC_PR_6837", "PBC_PR_6835_5A", "GBM"], ["A1_LN", "D1_LN"]),
    ("cross_holdout_Tonsil", ["A1_LN", "D1_LN", "Breast_Cancer", "PBC_PR_6837", "PBC_PR_6835_5A", "GBM"], ["A1_TNSL", "D1_TNSL", "HumanTonsil_10x", "HumanTonsil_AddOns_S1"]),
    ("cross_holdout_Breast_Cancer", ["A1_LN", "D1_LN", "A1_TNSL", "D1_TNSL", "HumanTonsil_10x", "HumanTonsil_AddOns_S1", "GBM"], ["Breast_Cancer", "PBC_PR_6837", "PBC_PR_6835_5A"]),
    ("cross_holdout_GBM", ["A1_LN", "D1_LN", "A1_TNSL", "D1_TNSL", "HumanTonsil_10x", "HumanTonsil_AddOns_S1", "Breast_Cancer", "PBC_PR_6837", "PBC_PR_6835_5A"], ["GBM"]),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_root", default=str(PROJECT_ROOT / "SCP-NET-work" / "runs" / "dpas_graph_10_dataset_20261101"))
    ap.add_argument("--epochs", type=int, default=100)
    ap.add_argument("--lr_scale", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--directions", default="all", help="all, or comma-separated direction names")
    args = ap.parse_args()
    genes, proteins = load_panel()
    selected = {x[0] for x in DIRECTIONS} if args.directions == "all" else {x.strip() for x in args.directions.split(",") if x.strip()}
    for name, sources, targets in DIRECTIONS:
        if name in selected:
            train_direction(name, sources, targets, args, genes, proteins)


if __name__ == "__main__":
    main()
