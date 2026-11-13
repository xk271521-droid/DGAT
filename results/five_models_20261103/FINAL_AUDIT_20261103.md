# 五模型比较表最终审计

更新日期：2026-11-03

## 覆盖与结构

- 五个模型各有 8 个同组织目标和 10 个跨组织目标，共 18 个目标/模型。
- 同组织页 40 行，跨组织页 50 行，总计 90 行。
- 每个目标数据集在对应实验页都有且仅有五个模型各一行；指标单元格及来源文件均完整，无重复记录或公式错误。
- 跨组织组固定为 LN、Tonsil、Breast Cancer、GBM；组内目标逐数据集单独评分。

## 修正版与数值核验

- ProteinConvex：表中使用新验证的 source-only RNA Neighbor Calibration。18/18 目标基线预测均以 `<=9.54e-7` 最大绝对误差重建；校准只用 source validation 选参数，不改 checkpoint、训练参数或数据。PCC-Pro、RMSE-Global、MAE-Global 在 18/18 目标均改善。未校准原始结果保留作消融。
- NicheTrans：Breast Cancer 跨组织 3 个目标来自严格 7-source 修正版；矩阵形状分别为 4160×31、4850×31、3243×31。
- DPAS-Graph：严格修正版 Breast Cancer 方向使用 7 个非乳腺癌 source，3 个 targets。prediction/truth 独立指标复算与 summary 的最大绝对差异为 `2.72e-15`。
- DGAT：12 个方向、18 个目标均核对 source/target 和数组形状。prediction/truth 独立复算与结果指标中的 PCC-Pro/RMSE-Global 最大绝对差异为 `1.68e-8`。
- ProteinConvex：使用已按模型内部蛋白顺序修正的 truth 矩阵；三个乳腺癌跨组织目标独立复算与结果记录一致。
- SPINE：使用 `analysis_ready` 最终结果；Breast Cancer 组来自 v6 正确整组留出结果。PBC 同组织出现接近零或负相关是 prediction/truth 独立复算一致的真实结果。

DPAS-Graph 修正版 Breast Cancer 方向 `best_epoch=90`、source validation PCC-Pro=`0.637846`。三个目标 PCC-Pro / RMSE-Global 为：Breast_Cancer `0.5446 / 0.3121`、PBC_PR_6837 `0.2531 / 0.2981`、PBC_PR_6835_5A `0.2661 / 0.3765`。

## 旧错误结果

NicheTrans 与 DPAS-Graph 旧 Breast Cancer source 泄漏结果已删除，不进入本表。当前 Breast Cancer 留出结果均用正确的 7-source/3-target 配置。

## 文件一致性

当前推荐工作簿 `five_model_comparison_20261103_rna_calibrated.xlsx` 的项目正式版、Codex 交付版和桌面版 SHA256 完全一致：`F410A1E26242E9F406A378A3B1D15747D34635AEEA08E8E623BA8CE51A9AB558`。原始未校准工作簿 `five_model_comparison_20261103.xlsx` 保留作基线，不覆盖。校准细节和对比见 `PROTEINCONVEX_RNA_CALIBRATION_10_DATASET_20261103.md`。
