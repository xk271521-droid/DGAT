# ProteinConvex 模型改进筛选协议

日期：2026-11-03

## 用途

降低每个 ProteinConvex 改进候选的初筛成本，同时保持训练数据隔离。这个两方向筛选分数用于候选筛选，不替代正式十数据集表中的 18 目标结果。

## 固定筛选方向

1. `cross_holdout_Breast_Cancer`：训练 source 为 `A1_LN`、`D1_LN`、`A1_TNSL`、`D1_TNSL`、`HumanTonsil_10x`、`HumanTonsil_AddOns_S1`、`GBM`；主评分 target 为 `PBC_PR_6837`。
2. `cross_holdout_GBM`：训练 source 为 `A1_LN`、`D1_LN`、`A1_TNSL`、`D1_TNSL`、`HumanTonsil_10x`、`HumanTonsil_AddOns_S1`、`Breast_Cancer`、`PBC_PR_6837`、`PBC_PR_6835_5A`；评分 target 为 `GBM`。

两项均为完整 source-to-target 训练与推理方向，不复用目标 ADT 选择 checkpoint 或校准参数。乳腺癌留出时不得把任何乳腺癌数据放入 source。

## 实验控制

- 保持十数据集、9384 个共同 RNA 基因、31 个共同 ADT 蛋白及其顺序不变。
- 每个候选只提出一个模型机制变化；seed、训练轮数、batch、optimizer、learning rate、loss 权重、source split 和预处理沿用正式 ProteinConvex 配置。若候选必须改这些数值，先单独记录为参数实验，不与结构候选混在一起。
- checkpoint 以及任何校准量只能由 source validation 选择；target ADT 只在最后生成评分时打开。
- 保存配置、模型代码快照、日志、checkpoint、target prediction/truth、逐蛋白指标和失败记录；每个候选用独立目录，不能覆盖正式基线。

## 评分和决策

首要筛选指标为两个 target 的 PCC-Pro 宏平均；PCC-Spot、Spearman-Pro/Spot 为相关性护栏，RMSE-Global、RMSE-Pro、MAE 为误差护栏。必须同时报告两个目标各自的指标，不得只报平均值。

任一目标有实质退化或误差护栏恶化，则记录为负结果并停止该候选。两项均稳定改善时，才扩展运行同组织 8 个方向及跨组织其余 targets，完整复核 18 个目标后才能更新正式比较表。

## 防止反复调测试集

重复根据同一两个 target 的指标改模型，会逐步把它们变成开发集。因此两方向结果只作候选筛选，不据此反复搜索参数。每个候选的设计和预期方向在运行前写入记录；进入论文主结果前，必须完成全目标评估，并把两方向筛选结果与剩余 target 分开报告。

## 当前校准版本参考基线

当前推荐基线为 `ProteinConvex + Strict RNA Neighbor Calibration`：

- `PBC_PR_6837`：PCC-Pro `0.247785`，PCC-Spot `0.230736`，RMSE-Global `0.270993`；
- `GBM`：PCC-Pro `0.204786`，PCC-Spot `0.205640`，RMSE-Global `0.313798`。

这些数值来自 source-only 校准版本；原始未校准 prediction 仍保留作消融对照。校准完整结果见 `SCP-NET-work/experiments/protein_convex_rna_calibration_10_dataset_20261103/RESULTS.md`。

每个新候选在 `SCP-NET-work/experiments/protein_convex_improvement_screening_20261103/` 下建立独立记录，结果放入对应独立 `runs` 子目录。
