#!/usr/bin/env bash
set -euo pipefail
source /extra_public/liqiang/miniconda3/etc/profile.d/conda.sh
conda activate HyperST
export PYTHONPATH=/public/liqiang/DGAT/upstream/SPINE
ROOT=/public/liqiang/DGAT/ST2P/data/processed/DGAT_final10_common9384_muon017_axis0_20261101
ADAPT=/public/liqiang/DGAT/SCP-NET-work/experiments/spine_10_dataset_20261101
SPINE=/public/liqiang/DGAT/upstream/SPINE
RUNROOT=/public/liqiang/DGAT/SCP-NET-work/runs/spine_10_dataset_20261101_v4
mkdir -p "$RUNROOT" "$RUNROOT/locks" "$ADAPT/logs"

run_one(){
  local gpu="$1" name="$2"; shift 2
  local lock="$RUNROOT/locks/$name.lock"
  if ! mkdir "$lock" 2>/dev/null; then echo "SKIP_LOCKED $name"; return 0; fi
  local lock="$RUNROOT/locks/$name.lock"
  if ! mkdir "$lock" 2>/dev/null; then echo "SKIP_LOCKED $name"; return 0; fi
  local sources=("$@")
  local targets=()
  case "$name" in
    same_LN_A1_to_D1) targets=(D1_LN);;
    same_LN_D1_to_A1) targets=(A1_LN);;
    same_TNSL_A1_to_D1) targets=(D1_TNSL);;
    same_TNSL_D1_to_A1) targets=(A1_TNSL);;
    same_HumanTonsil_10x_to_AddOns) targets=(HumanTonsil_AddOns_S1);;
    same_HumanTonsil_AddOns_to_10x) targets=(HumanTonsil_10x);;
    same_PBC_6837_to_6835) targets=(PBC_PR_6835_5A);;
    same_PBC_6835_to_6837) targets=(PBC_PR_6837);;
    cross_holdout_LN) targets=(A1_LN D1_LN);;
    cross_holdout_Tonsil) targets=(A1_TNSL D1_TNSL HumanTonsil_10x HumanTonsil_AddOns_S1);;
    cross_holdout_Breast_Cancer) targets=(Breast_Cancer PBC_PR_6837 PBC_PR_6835_5A);;
    cross_holdout_GBM) targets=(GBM);;
    *) echo "unknown job $name"; return 2;;
  esac
  local droot="$RUNROOT/data/$name"
  python "$ADAPT/prepare_spine_direction.py" --prepared_root "$ROOT" --spine_root "$droot" --name "$name" --sources "${sources[@]}" --targets "${targets[@]}" > "$RUNROOT/${name}.prepare.log" 2>&1
  CUDA_VISIBLE_DEVICES="$gpu" python "$ADAPT/train_rna_to_protein_adapted.py" \
    --dataset "$name" --source_dataroot "$droot/$name/dataset" --embed_dataroot "$droot/$name/dataset/embed_dataroot" \
    --save_dir "$RUNROOT/models" --protein_list protein_list.json --rna_gene_list rna_gene_list.json \
    --feature_dim 9384 --n_proteins 31 --batch_size 1 --gradient_accumulation_steps 2 --epochs 500 --sample_times 10 \
    --num_workers 1 --normalize_method None --rna_feature_normalize none --eval_step 1 --early_stop_patience 50 \
    --device 0 --exp_code "$name" --hidden_dim 256 --pairwise_hidden_dim 256 --n_layers 6 --dropout 0.3 \
    --attn_dropout 0.3 --n_heads 8 --mlp_intermediate_dim 1024 --mlp_num_layers 2 \
    > "$RUNROOT/${name}.train.log" 2>&1
  python "$ADAPT/evaluate_spine_targets.py" --spine_repo "$SPINE" --dataset_root "$droot/$name/dataset" \
    --embed_root "$droot/$name/dataset/embed_dataroot" --dataset "$name" --target_ids "${targets[@]}" \
    --model_dir "$RUNROOT/models/${name}/${name}" > "$RUNROOT/${name}.evaluate.log" 2>&1
}

run_one 2 same_LN_A1_to_D1 A1_LN &
run_one 3 same_LN_D1_to_A1 D1_LN &
wait
run_one 2 same_TNSL_A1_to_D1 A1_TNSL &
run_one 3 same_TNSL_D1_to_A1 D1_TNSL &
wait
run_one 2 same_HumanTonsil_10x_to_AddOns HumanTonsil_10x &
run_one 3 same_HumanTonsil_AddOns_to_10x HumanTonsil_AddOns_S1 &
wait
run_one 2 same_PBC_6837_to_6835 PBC_PR_6837 &
run_one 3 same_PBC_6835_to_6837 PBC_PR_6835_5A &
wait
run_one 2 cross_holdout_LN A1_TNSL D1_TNSL HumanTonsil_10x HumanTonsil_AddOns_S1 Breast_Cancer PBC_PR_6837 PBC_PR_6835_5A GBM &
run_one 3 cross_holdout_Tonsil A1_LN D1_LN Breast_Cancer PBC_PR_6837 PBC_PR_6835_5A GBM &
wait
run_one 2 cross_holdout_Breast_Cancer A1_LN D1_LN A1_TNSL D1_TNSL HumanTonsil_10x HumanTonsil_AddOns_S1 PBC_PR_6837 PBC_PR_6835_5A GBM &
run_one 3 cross_holdout_GBM A1_LN D1_LN A1_TNSL D1_TNSL HumanTonsil_10x HumanTonsil_AddOns_S1 Breast_Cancer PBC_PR_6837 PBC_PR_6835_5A &
wait
echo "SPINE_ALL_DONE"
