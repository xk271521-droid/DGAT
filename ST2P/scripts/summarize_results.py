#!/usr/bin/env python3
"""Create an auditable flat table from ST2P result directories."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


FIELDS = [
    "result_dir", "model", "target", "sources", "seed", "epochs", "best_epoch",
    "best_source_validation_pcc_pro_mean", "strict_pcc_pro", "strict_pcc_spot",
    "strict_spearman_pro", "strict_spearman_spot", "strict_rmse_global",
    "strict_rmse_pro", "strict_mae_global", "oracle_epoch", "oracle_pcc_pro",
    "registry", "selection",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    rows = []
    for path in sorted(args.results_root.rglob("results.json")):
        try:
            record = json.loads(path.read_text())
            strict = record.get("target_source_selected") or {}
        except (OSError, json.JSONDecodeError):
            continue
        rows.append({
            "result_dir": str(path.parent), "model": record.get("model"),
            "target": record.get("target"), "sources": ";".join(record.get("sources", [])),
            "seed": record.get("seed"), "epochs": record.get("epochs"),
            "best_epoch": record.get("best_epoch"),
            "best_source_validation_pcc_pro_mean": record.get("best_source_validation_pcc_pro_mean"),
            "strict_pcc_pro": strict.get("pcc_pro_mean"),
            "strict_pcc_spot": strict.get("pcc_spot_mean"),
            "strict_spearman_pro": strict.get("spearman_pro_mean"),
            "strict_spearman_spot": strict.get("spearman_spot_mean"),
            "strict_rmse_global": strict.get("rmse_global"),
            "strict_rmse_pro": record.get("target_source_selected_rmse_pro"),
            "strict_mae_global": strict.get("mae_global"),
            "oracle_epoch": record.get("best_test_oracle_epoch"),
            "oracle_pcc_pro": record.get("best_test_oracle_pcc_pro_mean"),
            "registry": record.get("preprocessing"), "selection": record.get("selection"),
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
