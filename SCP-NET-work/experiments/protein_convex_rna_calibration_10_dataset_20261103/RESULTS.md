# ProteinConvex 十数据集 RNA 邻居校准结果

日期：2026-11-03

## 方法边界

这是对十数据集 ProteinConvex/SCN-Pro 已训练 checkpoint 的严格 source-only RNA 邻居校准。模型结构、权重、训练参数、数据和 9384/31 面板均未改变，也没有重新训练。

每个正式方向只使用 source validation 的蛋白表达选择空间/RNA 图混合比和逐蛋白平滑强度。Target 图仅使用 RNA 与空间坐标；target ADT 在所有校准参数固定后才用于最终评分。每个方向先用保存的 checkpoint 重建原始预测，并与原保存 prediction 比较，18 个方向最大绝对差异为 `9.54e-7`。

## 结果

以下为 18 个目标等权宏平均。相关性越高越好，RMSE/MAE 越低越好。

| 实验范围 | 指标 | 原始 ProteinConvex | RNA 校准后 | 变化 |
|---|---|---:|---:|---:|
| 全部 18 目标 | PCC-Pro | 0.451864 | 0.463398 | +0.011534 |
| 全部 18 目标 | PCC-Spot | 0.483152 | 0.503344 | +0.020192 |
| 全部 18 目标 | Spearman-Pro | 0.409382 | 0.421741 | +0.012359 |
| 全部 18 目标 | Spearman-Spot | 0.451718 | 0.463352 | +0.011635 |
| 全部 18 目标 | RMSE-Global | 0.323224 | 0.307327 | -0.015897 |
| 全部 18 目标 | MAE-Global | 0.205049 | 0.197128 | -0.007920 |
| 同组织 8 目标 | PCC-Pro | 0.517581 | 0.534085 | +0.016504 |
| 同组织 8 目标 | PCC-Spot | 0.571922 | 0.595889 | +0.023968 |
| 同组织 8 目标 | RMSE-Global | 0.291529 | 0.273489 | -0.018040 |
| 跨组织 10 目标 | PCC-Pro | 0.399290 | 0.406848 | +0.007558 |
| 跨组织 10 目标 | PCC-Spot | 0.412136 | 0.429308 | +0.017172 |
| 跨组织 10 目标 | RMSE-Global | 0.348580 | 0.334397 | -0.014183 |

逐目标对比中，PCC-Pro、RMSE-Global、MAE-Global 均为 18/18 目标改善；PCC-Spot 为 16/18、Spearman-Pro 为 18/18、Spearman-Spot 为 17/18 改善。所有方向使用同一预先固定的 `rho`、`lambda` 网格，没有依据 target ADT 改网格或挑选方向。

## 结论与限制

该候选在当前十数据集协议的全部 18 个目标上，对核心 Pearson-Pro 与误差指标呈一致改善，值得作为 ProteinConvex 的严格 source-only 校准版本继续使用。它是推理阶段的 RNA/空间邻居校准，不是重训后的新网络；论文和表格应与未校准 ProteinConvex 基线区分并报告两者。

校准参数与 checkpoint 都使用同一方向的 source validation，因此虽然没有 target ADT 泄漏，source validation 被重复用于模型/校准选择仍有过拟合风险。正式论文若需要更严格的泛化主张，后续可用 source 组织内层留一或嵌套验证再评估；本次校准测试集没有参与参数选择。

逐方向参数、source-validation 网格分数、原始/校准指标和校准 prediction 位于 `SCP-NET-work/runs/protein_convex_rna_calibration_10_dataset_20261103/all_directions/`。原始基线 prediction 和 checkpoint 保留在 `SCP-NET-work/runs/protein_convex_10_dataset_20261101/`。
