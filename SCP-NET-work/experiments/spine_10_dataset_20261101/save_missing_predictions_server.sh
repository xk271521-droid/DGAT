#!/usr/bin/env bash
set -euo pipefail
source /extra_public/liqiang/miniconda3/etc/profile.d/conda.sh
conda activate HyperST
export PYTHONPATH=/public/liqiang/DGAT/upstream/SPINE
ADAPT=/public/liqiang/DGAT/SCP-NET-work/experiments/spine_10_dataset_20261101
SPINE=/public/liqiang/DGAT/upstream/SPINE
RUNROOT=/public/liqiang/DGAT/SCP-NET-work/runs/spine_10_dataset_20261101_v5
SCRIPT=/public/liqiang/DGAT/SCP-NET-work/experiments/spine_10_dataset_20261101/run_remaining_server_v5.sh

# Let the queued training/evaluation sequence finish before filling any gaps.
while pgrep -f '[r]un_remaining_server_v5.sh' >/dev/null; do sleep 60; done

evaluate(){
  local name="$1"; shift
  local targets=("$@")
  local droot="$RUNROOT/data/$name/$name"
  local model="$RUNROOT/models/$name/$name"
  if [ ! -f "$model/split0/checkpoints/best.pth" ]; then
    echo "MISSING_CHECKPOINT $name"; return 1
  fi
  CUDA_VISIBLE_DEVICES='' python "$ADAPT/evaluate_spine_targets.py" --gpu -1 \
    --spine_repo "$SPINE" --dataset_root "$droot/dataset" \
    --embed_root "$droot/dataset/embed_dataroot" --dataset "$name" \
    --target_ids "${targets[@]}" --model_dir "$model" \
    > "$RUNROOT/${name}.evaluate_cpu.log" 2>&1
}

evaluate same_HumanTonsil_10x_to_AddOns HumanTonsil_AddOns_S1
evaluate same_HumanTonsil_AddOns_to_10x HumanTonsil_10x
evaluate cross_holdout_Breast_Cancer Breast_Cancer PBC_PR_6837 PBC_PR_6835_5A
evaluate cross_holdout_GBM GBM
echo "ALL_HELDOUT_PREDICTIONS_SAVED"
