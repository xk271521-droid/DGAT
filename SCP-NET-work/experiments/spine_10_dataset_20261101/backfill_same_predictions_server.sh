#!/usr/bin/env bash
set -euo pipefail
source /extra_public/liqiang/miniconda3/etc/profile.d/conda.sh
conda activate HyperST
export PYTHONPATH=/public/liqiang/DGAT/upstream/SPINE
ADAPT=/public/liqiang/DGAT/SCP-NET-work/experiments/spine_10_dataset_20261101
SPINE=/public/liqiang/DGAT/upstream/SPINE
BASE=/public/liqiang/DGAT/SCP-NET-work/runs

eval_one(){
  local runroot="$1" name="$2"; shift 2
  local targets=("$@")
  local droot="$runroot/data/$name/$name"
  local model="$runroot/models/$name/$name"
  local out="$model/heldout_predictions"
  if [ -f "$out/$(printf %s "${targets[0]}")_predictions.h5ad" ]; then
    echo "SKIP_ALREADY_SAVED $name"; return 0
  fi
  test -f "$model/split0/checkpoints/best.pth" || { echo "MISSING_CHECKPOINT $name"; return 1; }
  CUDA_VISIBLE_DEVICES='' python "$ADAPT/evaluate_spine_targets.py" --gpu -1 \
    --spine_repo "$SPINE" --dataset_root "$droot/dataset" \
    --embed_root "$droot/dataset/embed_dataroot" --dataset "$name" \
    --target_ids "${targets[@]}" --model_dir "$model" \
    > "$runroot/${name}.backfill.log" 2>&1
  echo "SAVED $name"
}

eval_one "$BASE/spine_10_dataset_20261101_v3" same_LN_A1_to_D1 D1_LN
eval_one "$BASE/spine_10_dataset_20261101_v2" same_TNSL_D1_to_A1 A1_TNSL
eval_one "$BASE/spine_10_dataset_20261101_v2" same_PBC_6837_to_6835 PBC_PR_6835_5A
eval_one "$BASE/spine_10_dataset_20261101_v2" same_PBC_6835_to_6837 PBC_PR_6837
echo "SAME_ORGANIZATION_PREDICTIONS_COMPLETE"
