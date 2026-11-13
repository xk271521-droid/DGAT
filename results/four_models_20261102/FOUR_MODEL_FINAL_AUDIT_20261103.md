# 四模型表最终复核

复核日期：2026-11-03

## 当前纳入范围

本表只比较 ProteinConvex + Strict RNA、NicheTrans、SPINE、DGAT。DPAS-Graph 按用户要求暂不纳入；其旧 Breast Cancer 泄漏结果已删除。四个纳入模型各有 8 个同组织目标和 10 个跨组织目标，总计 72 条指标明细。跨组织仍按 LN、Tonsil、Breast Cancer、GBM 四组，每个目标数据集分别评分。

## DGAT 核验

- 正式结果目录：`SCP-NET-work/runs/dgat_10_dataset_20261101_all_results`
- 方向和目标覆盖：12 个训练方向、18 个 prediction/truth 目标文件
- 所有预测/真实值矩阵均为对应 spot 数 × 31 个蛋白，无非有限值
- source/target 划分逐方向核对通过；Breast Cancer 留出 source 为 `A1_LN`、`D1_LN`、`A1_TNSL`、`D1_TNSL`、`HumanTonsil_10x`、`HumanTonsil_AddOns_S1`、`GBM`
- 从 prediction/truth 独立重算 PCC-Pro、PCC-Spot、Spearman-Pro、Spearman-Spot、RMSE、MAE；与保存指标中的 PCC-Pro 和 RMSE-Global 最大绝对差异为 `1.68e-8`

## NicheTrans 修正版

乳腺癌留出三条结果来自 `SCP-NET-work/runs/nichetrans_10_dataset_20261102_corrected_breast_holdout/cross_holdout_Breast_Cancer`，source 为上述 7 个非乳腺癌数据集，targets 为 `Breast_Cancer`、`PBC_PR_6837`、`PBC_PR_6835_5A`。三组矩阵形状分别为 4160×31、4850×31、3243×31。指标从修正版 prediction/truth 重新计算后写入表格。

## 文件与排除项

结果表为 `comparison/four_model_comparison_20261102.xlsx`。项目结果副本、交付副本和桌面副本的 SHA256 经核对一致。工作表数据行中不含 DPAS-Graph 旧值；DPAS-Graph 等待严格 Breast Cancer 修正版，不影响当前四模型表的 72 条覆盖。
