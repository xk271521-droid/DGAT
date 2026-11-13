#!/usr/bin/env python3
"""NicheTrans official model on the frozen ten-dataset registry.

The official model source is imported unchanged from upstream/NicheTrans.
Only the data adapter, common-panel dimensions, source-only validation, and
experiment grouping are implemented here.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.stats
import torch
from sklearn.neighbors import NearestNeighbors
from torch.utils.data import ConcatDataset, DataLoader, Dataset, Subset

UPSTREAM = Path(r"E:\DGAT\upstream\NicheTrans")
PANEL_PATH = Path(r"E:\DGAT\SCP-NET-work\docs\ten_dataset_common_panel_20261101.json")
sys.path.insert(0, str(UPSTREAM))
from model.nicheTrans import NicheTrans  # noqa: E402


class NicheDataset(Dataset):
    def __init__(self, x: np.ndarray, y: np.ndarray, neighbors: np.ndarray):
        self.x = torch.from_numpy(x.astype(np.float32, copy=False))
        self.y = torch.from_numpy(y.astype(np.float32, copy=False))
        self.neighbors = torch.from_numpy(neighbors.astype(np.float32, copy=False))

    def __len__(self):
        return self.x.shape[0]

    def __getitem__(self, i):
        return self.x[i], self.y[i], self.neighbors[i]


def metrics(pred: np.ndarray, truth: np.ndarray) -> dict:
    pcc, spearman, rmse = [], [], []
    for j in range(truth.shape[1]):
        if np.std(truth[:, j]) <= 1e-8 or np.std(pred[:, j]) <= 1e-8:
            continue
        pcc.append(float(scipy.stats.pearsonr(pred[:, j], truth[:, j]).statistic))
        spearman.append(float(scipy.stats.spearmanr(pred[:, j], truth[:, j]).statistic))
        rmse.append(float(np.sqrt(np.mean((pred[:, j] - truth[:, j]) ** 2))))
    return {
        "pcc_pro": float(np.mean(pcc)) if pcc else float("nan"),
        "spearman_pro": float(np.mean(spearman)) if spearman else float("nan"),
        "rmse_global": float(np.sqrt(np.mean((pred - truth) ** 2))),
        "rmse_pro": float(np.mean(rmse)) if rmse else float("nan"),
        "n_valid_proteins": len(pcc),
    }


def load_section(root: Path, name: str, genes: list[str], proteins: list[str], k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rna = ad.read_h5ad(root / name / "rna_dgat.h5ad")[:, genes].copy()
    adt = ad.read_h5ad(root / name / "adt_dgat_clr.h5ad")[:, proteins].copy()
    if not rna.obs_names.equals(adt.obs_names):
        raise ValueError(f"RNA/ADT barcode mismatch: {name}")
    coords = np.asarray(rna.obsm["spatial"], dtype=np.float32)
    x = np.asarray(rna.X.toarray() if hasattr(rna.X, "toarray") else rna.X, dtype=np.float32)
    y = np.asarray(adt.X.toarray() if hasattr(adt.X, "toarray") else adt.X, dtype=np.float32)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords)
    idx = idx[:, 1:]
    neighbors = x[idx]
    return x, y, neighbors


def predict(model, loader, device):
    model.eval()
    pred, truth = [], []
    with torch.no_grad():
        for x, y, n in loader:
            pred.append(model(x.to(device), n.to(device)).cpu().numpy())
            truth.append(y.numpy())
    return np.concatenate(pred), np.concatenate(truth)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepared_root", required=True, type=Path)
    ap.add_argument("--out_root", required=True, type=Path)
    ap.add_argument("--name", required=True)
    ap.add_argument("--sources", nargs="+", required=True)
    ap.add_argument("--targets", nargs="+", required=True)
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--batch_size", type=int, default=32)
    ap.add_argument("--neighbors", type=int, default=12)
    ap.add_argument("--noise_rate", type=float, default=0.2)
    ap.add_argument("--dropout_rate", type=float, default=0.2)
    ap.add_argument("--lr", type=float, default=0.0003)
    ap.add_argument("--weight_decay", type=float, default=5e-4)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    random.seed(args.seed); np.random.seed(args.seed); torch.manual_seed(args.seed)
    if torch.cuda.is_available(): torch.cuda.manual_seed_all(args.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    out = args.out_root / args.name
    if out.exists():
        raise FileExistsError(out)
    out.mkdir(parents=True)

    manifest = json.loads((args.prepared_root / "prepared_manifest.json").read_text(encoding="utf-8"))
    names = [x["dataset"] for x in manifest["datasets"]]
    panel = json.loads(PANEL_PATH.read_text(encoding="utf-8"))
    genes = [str(x) for x in panel["rna_common_genes"]]
    proteins = [str(x) for x in panel["adt_common_proteins"]]
    if len(genes) != 9384 or len(proteins) != 31:
        raise ValueError(f"Unexpected ten-dataset panel: genes={len(genes)}, proteins={len(proteins)}")

    blocks = {}
    for n in [*args.sources, *args.targets]:
        blocks[n] = load_section(args.prepared_root, n, genes, proteins, args.neighbors)

    datasets, fit_indices, val_indices = [], [], []
    for pos, n in enumerate(args.sources):
        x, y, neigh = blocks[n]
        ds = NicheDataset(x, y, neigh)
        datasets.append(ds)
        order = np.arange(len(ds)); np.random.default_rng(args.seed + pos).shuffle(order)
        n_val = max(1, int(round(len(order) * 0.1)))
        val_indices.append((len(datasets) - 1, order[:n_val]))
        fit_indices.append((len(datasets) - 1, order[n_val:]))
    train_ds = ConcatDataset(datasets)
    def global_indices(parts):
        offsets=[]; cur=0
        for d in datasets: offsets.append(cur); cur += len(d)
        return np.concatenate([offsets[i] + idx for i, idx in parts]).astype(np.int64)
    train_loader = DataLoader(Subset(train_ds, global_indices(fit_indices).tolist()), batch_size=args.batch_size, shuffle=True, drop_last=True, num_workers=0)
    val_loader = DataLoader(Subset(train_ds, global_indices(val_indices).tolist()), batch_size=args.batch_size, shuffle=False, num_workers=0)

    model = NicheTrans(source_length=len(genes), target_length=len(proteins), noise_rate=args.noise_rate, dropout_rate=args.dropout_rate).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    criterion = torch.nn.MSELoss()
    best = -float("inf"); best_epoch = 0; history=[]
    for epoch in range(1, args.epochs + 1):
        model.train(); losses=[]
        for x, y, neigh in train_loader:
            x, y, neigh = x.to(device), y.to(device), neigh.to(device)
            if random.random() > 0.7:
                neigh = neigh * torch.bernoulli(torch.full((neigh.size(0), neigh.size(1), 1), 0.5, device=device))
            pred = model(x, neigh); loss = criterion(pred, y)
            optimizer.zero_grad(); loss.backward(); optimizer.step(); losses.append(float(loss.detach().cpu()))
        vp, vt = predict(model, val_loader, device); vm = metrics(vp, vt)
        row={"epoch":epoch,"train_loss":float(np.mean(losses)),**{f"val_{k}":v for k,v in vm.items()}}; history.append(row)
        print(f"[{args.name}] epoch={epoch:02d}/{args.epochs} loss={row['train_loss']:.5f} val_pcc={vm['pcc_pro']:.4f}",flush=True)
        if vm["pcc_pro"] > best:
            best=vm["pcc_pro"]; best_epoch=epoch; torch.save(model.state_dict(),out/"best_source_validation.pth")
    (out/"metrics_history.json").write_text(json.dumps(history,indent=2),encoding="utf-8")
    model.load_state_dict(torch.load(out/"best_source_validation.pth",map_location=device,weights_only=True))
    summary={"name":args.name,"sources":args.sources,"targets":args.targets,"genes":len(genes),"proteins":proteins,"best_epoch":best_epoch,"best_source_validation_pcc_pro":best,"official_parameters":{"noise_rate":args.noise_rate,"dropout_rate":args.dropout_rate,"epochs":args.epochs,"batch_size":args.batch_size,"neighbors":args.neighbors,"lr":args.lr,"weight_decay":args.weight_decay,"seed":args.seed},"device":str(device),"target_metrics":{}}
    for n in args.targets:
        x,y,neigh=blocks[n]; pred,truth=predict(model,DataLoader(NicheDataset(x,y,neigh),batch_size=args.batch_size,shuffle=False),device); m=metrics(pred,truth); summary["target_metrics"][n]=m; np.save(out/f"{n}_predictions.npy",pred); np.save(out/f"{n}_truth.npy",truth); (out/f"{n}_metrics.json").write_text(json.dumps(m,indent=2),encoding="utf-8")
    (out/"results.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(f"[{args.name}] COMPLETE best_epoch={best_epoch} targets={list(summary['target_metrics'])}",flush=True)


if __name__ == "__main__":
    main()
