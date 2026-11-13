from __future__ import annotations

import argparse
from argparse import Namespace
from pathlib import Path

from run_dgat_10_dataset import load_panel, run_one


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepared_root", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out_root", required=True)
    ap.add_argument("--epochs", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    sources = [
        "A1_LN",
        "D1_LN",
        "A1_TNSL",
        "D1_TNSL",
        "HumanTonsil_10x",
        "HumanTonsil_AddOns_S1",
        "GBM",
    ]
    targets = ["Breast_Cancer", "PBC_PR_6837", "PBC_PR_6835_5A"]
    genes, proteins = load_panel(args.panel)
    run_args = Namespace(
        prepared_root=args.prepared_root,
        out_root=args.out_root,
        epochs=args.epochs,
        seed=args.seed,
    )
    run_one(
        "cross_holdout_Breast_Cancer",
        sources,
        targets,
        run_args,
        genes,
        proteins,
    )


if __name__ == "__main__":
    main()
