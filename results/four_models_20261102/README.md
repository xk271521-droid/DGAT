# 四模型十数据集比较表（已被五模型版取代）

整理日期：2026-11-02

当前比较模型：ProteinConvex + Strict RNA、NicheTrans、SPINE、DGAT。DPAS-Graph 的 Breast Cancer 留出方向尚待严格重跑，因此本版不纳入 DPAS-Graph。

本目录中的四模型工作簿是中间版本，已于 2026-11-03 被 `results/five_models_20261103/comparison/five_model_comparison_20261103.xlsx` 取代。后续使用五模型目录中的正式表。

统一输入：10 个数据集、9384 个共同 RNA 基因、31 个共同 ADT 蛋白。

实验组织：

- 同组织：LN、TNSL、Human Tonsil、Breast Cancer 四组配对数据，8 个双向方向；
- 跨组织：留出 LN、Tonsil、Breast Cancer、GBM 四组。每组下面按目标数据集单独评分，不能把同一留出组下的多个目标合并成额外组别。

正式比较表位于 `comparison/four_model_comparison_20261102.xlsx`。每个明细行保存目标数据集、来源方向、模型、8 个指标、Best epoch、结果状态和来源路径。

SPINE 使用 `results/SPINE_10_dataset_20261101/analysis_ready` 中的最终 prediction/truth，按相同指标重新计算。DGAT 已从 `dgat_10_dataset_20261101_all_results` 拉回，覆盖 12 个方向和 18 个目标；其中 Breast Cancer 方向使用 GPU2 修正版。

NicheTrans Breast Cancer 留出已使用正确的 7-source/3-target 修正版。最终表包含 4 个模型、72 条有效明细；DPAS-Graph 的旧 Breast Cancer 泄漏结果已删除，按用户要求暂不纳入本版。DGAT 的 18 个目标结果已逐项核对协议、预测/真实值和指标。协议与数值复核见 `FOUR_MODEL_FINAL_AUDIT_20261103.md`。
