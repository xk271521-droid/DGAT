#!/usr/bin/env python3
"""Audit source-only RNA calibration arrays and metrics against frozen runs."""
import json
from pathlib import Path

import numpy as np
import sys

BASE = Path(r"E:/DGAT")
RUNS = BASE / "SCP-NET-work/runs"
OUT = RUNS / "protein_convex_rna_calibration_10_dataset_20261103/all_directions"
sys.path.insert(0, str(BASE / "ST2P/src"))
from triscale_center_protein_v9_dgat_compatible import dgat_metrics

summary = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))
if summary["direction_count"] != 18 or len(summary["records"]) != 18:
    raise RuntimeError("calibration must cover all 18 formal targets")

errors = []
same_count = cross_count = 0
improved = {k: 0 for k in ["pcc_pro_mean", "pcc_spot_mean", "spearman_pro_mean", "spearman_spot_mean", "rmse_global", "mae_global"]}
sum_base, sum_cal = {}, {}
max_reconstruction_diff = 0.0
for record in summary["records"]:
    direction = record["direction"]
    target = record["target"]
    base_dir = RUNS / "protein_convex_10_dataset_20261101" / direction
    candidate_dir = OUT / direction
    result = json.loads((base_dir / "results.json").read_text(encoding="utf-8"))
    if record["sources"] != result["sources"] or record["target"] != result["target"]:
        raise RuntimeError(f"source/target mismatch: {direction}")
    if not record["source_validation_only_selection"] or not record["target_adt_opened_after_selection"]:
        raise RuntimeError(f"selection protocol flags missing: {direction}")
    p = np.load(candidate_dir / "rna_calibrated_predictions.npy")
    y = np.load(base_dir / "target_source_selected_truth.npy")
    baseline_p = np.load(base_dir / "target_source_selected_predictions.npy")
    if p.shape != y.shape or p.shape != baseline_p.shape or p.shape[1] != 31:
        raise RuntimeError(f"array shape mismatch: {direction} {p.shape}/{y.shape}/{baseline_p.shape}")
    if not np.isfinite(p).all() or not np.isfinite(y).all():
        raise RuntimeError(f"non-finite array: {direction}")
    max_reconstruction_diff = max(max_reconstruction_diff, record["baseline_reconstruction_max_abs_diff"])
    if record["baseline_reconstruction_max_abs_diff"] > 2e-5:
        raise RuntimeError(f"raw prediction reproduction failed: {direction}")
    recomputed = dgat_metrics(p, y, record["common_panel"]["proteins"])
    for key in ["pcc_pro_mean", "pcc_spot_mean", "spearman_pro_mean", "spearman_spot_mean", "rmse_global", "mae_global"]:
        stored = record["calibrated_metrics"][key]
        if abs(float(recomputed[key]) - float(stored)) > 2e-5:
            raise RuntimeError(f"calibrated metric mismatch: {direction} {key}")
        sum_cal[key] = sum_cal.get(key, 0.0) + float(stored)
        base_value = float(record["baseline_metrics"][key])
        sum_base[key] = sum_base.get(key, 0.0) + base_value
        if key in ["rmse_global", "mae_global"]:
            improved[key] += int(float(stored) < base_value)
        else:
            improved[key] += int(float(stored) > base_value)
    if direction.startswith("same_"): same_count += 1
    else: cross_count += 1

if same_count != 8 or cross_count != 10:
    raise RuntimeError(f"wrong protocol counts: same={same_count}, cross={cross_count}")

report = {
    "direction_count": len(summary["records"]),
    "same_tissue_targets": same_count,
    "cross_tissue_targets": cross_count,
    "max_baseline_reconstruction_abs_diff": max_reconstruction_diff,
    "target_metric_improvements": {k: f"{v}/18" for k, v in improved.items()},
    "macro_mean": {k: {"baseline": sum_base[k] / 18, "calibrated": sum_cal[k] / 18, "delta": (sum_cal[k] - sum_base[k]) / 18} for k in sum_base},
    "status": "PASS",
}
(OUT / "AUDIT.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps(report, indent=2))
