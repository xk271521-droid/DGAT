# 五模型十数据集比较

整理日期：2026-11-03

当前纳入模型：ProteinConvex + source-only Strict RNA Neighbor Calibration、NicheTrans、DPAS-Graph、SPINE、DGAT。未经校准的 ProteinConvex 原始预测保留作基线；本表 ProteinConvex 行采用校准后预测。

统一实验协议：10 个数据集、9384 个共同 RNA 基因、31 个共同 ADT 蛋白；同组织 8 个双向目标，跨组织 4 个整组留出组、10 个目标。每个模型 18 个目标结果，共 90 条明细。

跨组织按 LN、Tonsil、Breast Cancer、GBM 分组；每组下每个测试数据集单独列出五个模型的指标。Breast Cancer 留出方向使用 NicheTrans、DPAS-Graph、DGAT 的 7-source/3-target 正确结果；不包含旧 PBC source 泄漏结果。

推荐工作簿：`comparison/five_model_comparison_20261103_rna_calibrated.xlsx`。原始未校准五模型表保留为 `comparison/five_model_comparison_20261103.xlsx`，用于基线/消融对照。完整审计记录见 `FINAL_AUDIT_20261103.md`。

ProteinConvex 校准逐方向参数、原始/校准指标和 prediction 位于 `SCP-NET-work/runs/protein_convex_rna_calibration_10_dataset_20261103/all_directions/summary.json` 及其方向子目录。校准实验汇总见 `PROTEINCONVEX_RNA_CALIBRATION_10_DATASET_20261103.md`。
