#!/usr/bin/env bash
set -euo pipefail
source /public/liqiang/miniconda3/etc/profile.d/conda.sh
conda activate HyperST
export CUDA_VISIBLE_DEVICES=4
export PYTHONUNBUFFERED=1
ROOT=/public/liqiang/DGAT
SCRIPT=$ROOT/SCP-NET-work/experiments/dgat_10_dataset_20261101/run_dgat_10_dataset.py
DATA=$ROOT/ST2P/data/processed/DGAT_final10_common9384_muon017_axis0_20261101
PANEL=$ROOT/SCP-NET-work/docs/ten_dataset_common_panel_20261101.json
OUT=$ROOT/SCP-NET-work/runs/dgat_10_dataset_20261101_server
mkdir -p "$OUT"
DIRECTIONS='all'
echo "START $(date -Is) CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
exec python "$SCRIPT" --prepared_root "$DATA" --panel "$PANEL" --out_root "$OUT" --epochs 100 --seed 0 --directions "$DIRECTIONS"
