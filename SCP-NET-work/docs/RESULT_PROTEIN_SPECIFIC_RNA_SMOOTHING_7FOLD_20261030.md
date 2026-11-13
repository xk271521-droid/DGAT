# 蛋白特异 RNA 相似邻居推理校准：早期探索结果（已被严格版本替代）

> 注意：本文件中的参数曾用其他外层折结果做快速探索，存在间接使用当前目标切片训练信息的风险。不要把本文件数字作为论文最终结果；请使用 `RESULT_PROTEIN_SPECIFIC_RNA_SMOOTHING_STRICT_7FOLD_20261030.md`。

日期：2026-10-30。

该方法不修改 ProteinConvex 训练和 checkpoint。对每个目标折：

1. 使用保存的 ProteinConvex 目标预测；
2. 目标切片只用 RNA 构建 8-NN 的 RNA 相似邻居图；
3. 每个蛋白单独选择平滑强度 `lambda_p`；
4. `lambda_p` 和图混合策略只用其他六个外层目标折的标签选择；
5. 当前目标 ADT 只用于最终评估，不参与参数选择。

交叉选择结果偏向 RNA 相似图（`rho=0`），而不是空间距离图。预测校准形式为：

```text
p'_i,p = (1 - lambda_p) * p_i,p
          + lambda_p * sum_j W_RNA(i,j) * p_j,p
```

## 七折结果

| 指标 | ProteinConvex 原预测 | 蛋白特异 RNA 校准 |
|---|---:|---:|
| PCC-Pro | 0.586593 | 0.594368 |
| PCC-Spot | 0.629418 | 0.645400 |
| RMSE-Global | 0.285059 | 0.272326 |

## Breast_Cancer 折

| 指标 | ProteinConvex 原预测 | 蛋白特异 RNA 校准 |
|---|---:|---:|
| PCC-Pro | 0.602035 | 0.611970 |
| PCC-Spot | 0.633267 | 0.642120 |
| RMSE-Global | 0.299574 | 0.291290 |

该结果是当前项目中第一个同时改善七折均值和 Breast_Cancer 单折三项主要指标的新增方案。它属于推理期校准，不改变 ProteinConvex 模型参数；正式报告应同时列出原预测和校准后预测。

复现脚本：`SCP-NET-work/experiments/spatial_hybrid_smoothing_20261030/analyze_spatial_hybrid_smoothing.py`。
逐折结果：`SCP-NET-work/experiments/spatial_hybrid_smoothing_20261030/results.csv`。
