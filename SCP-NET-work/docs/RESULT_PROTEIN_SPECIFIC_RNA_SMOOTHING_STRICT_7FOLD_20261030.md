# 蛋白特异 RNA 相似邻居推理校准：严格源验证七折结果

这是对前一版探索性结果的严格重跑。每个外层目标切片的校准参数只使用该折六个源切片的 validation spots 选择；当前目标切片只使用 RNA 建图和预测，目标 ADT 只在最终评估阶段读取。

对每个蛋白单独选择 `lambda_p`，并在空间邻居与 RNA 相似邻居之间选择图权重。最终选择均为 RNA 相似图主导（`rho=0`）。

## 七折结果

| 指标 | ProteinConvex 原预测 | 严格 RNA 校准 |
|---|---:|---:|
| PCC-Pro | 0.586593 | 0.594800 |
| PCC-Spot | 0.629418 | 0.644555 |
| RMSE-Global | 0.285059 | 0.272659 |

## Breast_Cancer 折

| 指标 | ProteinConvex 原预测 | 严格 RNA 校准 |
|---|---:|---:|
| PCC-Pro | 0.602035 | 0.612301 |
| PCC-Spot | 0.633267 | 0.642322 |
| RMSE-Global | 0.299574 | 0.290963 |

该结果遵守当前七切片 LODO 的 source-only 选择要求，可以作为完整流程结果报告。它仍属于 ProteinConvex 的推理校准版本，论文中应与原始 ProteinConvex 做消融对照。

严格复现脚本：`SCP-NET-work/experiments/strict_spatial_calibration_20261030/strict_spatial_calibration.py`。
逐折结果：`SCP-NET-work/experiments/strict_spatial_calibration_20261030/results_strict_source_validation.csv`。
