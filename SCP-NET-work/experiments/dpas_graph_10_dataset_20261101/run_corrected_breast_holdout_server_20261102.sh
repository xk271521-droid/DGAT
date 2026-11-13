#!/usr/bin/env bash
set -euo pipefail
export CUDA_VISIBLE_DEVICES=0
PY=/extra_public/liqiang/miniconda3/envs/HyperST/bin/python
SCRIPT=/public/liqiang/DGAT/server_runs/dpas_graph_20261101/dpas_graph_10_dataset_20261101/run_dpas_graph_10_dataset.py
OUT=/public/liqiang/DGAT/SCP-NET-work/runs/dpas_graph_10_dataset_20261102_corrected_breast_holdout
mkdir -p "$OUT"
exec "$PY" "$SCRIPT" --out_root "$OUT" --epochs 100 --seed 0 --directions cross_holdout_Breast_Cancer
