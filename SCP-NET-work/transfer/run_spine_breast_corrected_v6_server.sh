#!/usr/bin/env bash
set -euo pipefail

source /extra_public/liqiang/miniconda3/etc/profile.d/conda.sh
conda activate HyperST
export PYTHONPATH=/public/liqiang/DGAT/upstream/SPINE

ROOT=/public/liqiang/DGAT/ST2P/data/processed/DGAT_final10_common9384_muon017_axis0_20261101
ADAPT=/public/liqiang/DGAT/SCP-NET-work/experiments/spine_10_dataset_20261101
SPINE=/public/liqiang/DGAT/upstream/SPINE
RUNROOT=/public/liqiang/DGAT/SCP-NET-work/runs/spine_10_dataset_20261101_v6
NAME=cross_holdout_Breast_Cancer
GPU=3

SOURCES=(A1_LN D1_LN A1_TNSL D1_TNSL HumanTonsil_10x HumanTonsil_AddOns_S1 GBM)
TARGETS=(Breast_Cancer PBC_PR_6837 PBC_PR_6835_5A)

mkdir -p "$RUNROOT"
python "$ADAPT/prepare_spine_direction.py" \
  --prepared_root "$ROOT" \
  --spine_root "$RUNROOT/data/$NAME" \
  --name "$NAME" \
  --sources "${SOURCES[@]}" \
  --targets "${TARGETS[@]}" \
  > "$RUNROOT/$NAME.prepare.log" 2>&1

CUDA_VISIBLE_DEVICES="$GPU" python "$ADAPT/train_rna_to_protein_adapted.py" \
  --dataset "$NAME" \
  --source_dataroot "$RUNROOT/data/$NAME/$NAME/dataset" \
  --embed_dataroot "$RUNROOT/data/$NAME/$NAME/dataset/embed_dataroot" \
  --save_dir "$RUNROOT/models" \
  --protein_list protein_list.json \
  --rna_gene_list rna_gene_list.json \
  --feature_dim 9384 \
  --n_proteins 31 \
  --batch_size 1 \
  --gradient_accumulation_steps 2 \
  --epochs 500 \
  --sample_times 10 \
  --num_workers 1 \
  --normalize_method None \
  --rna_feature_normalize none \
  --eval_step 1 \
  --early_stop_patience 50 \
  --device 0 \
  --exp_code "$NAME" \
  --hidden_dim 256 \
  --pairwise_hidden_dim 256 \
  --n_layers 6 \
  --dropout 0.3 \
  --attn_dropout 0.3 \
  --n_heads 8 \
  --mlp_intermediate_dim 1024 \
  --mlp_num_layers 2 \
  > "$RUNROOT/$NAME.train.log" 2>&1

CUDA_VISIBLE_DEVICES="$GPU" python "$ADAPT/evaluate_spine_targets.py" \
  --gpu 0 \
  --spine_repo "$SPINE" \
  --dataset_root "$RUNROOT/data/$NAME/$NAME/dataset" \
  --embed_root "$RUNROOT/data/$NAME/$NAME/dataset/embed_dataroot" \
  --dataset "$NAME" \
  --target_ids "${TARGETS[@]}" \
  --model_dir "$RUNROOT/models/$NAME/$NAME" \
  > "$RUNROOT/$NAME.evaluate.log" 2>&1

echo "SPINE_CORRECTED_BREAST_CANCER_COMPLETE"
