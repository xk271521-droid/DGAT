# ProteinConvex 静态尺度 simplex 约束：Breast_Cancer

只把 CognateAnchor 的静态三尺度权重改为 softmax simplex，初始化仍为 `[0.5, 0.3, 0.2]`；其余模型、数据、seed、20 epoch、损失和 source-only 选点保持不变。

| 指标 | ProteinConvex 基线 | simplex 约束 |
|---|---:|---:|
| 源验证最佳 PCC-Pro | 0.751729 | 0.750758 |
| Breast PCC-Pro | 0.602035 | 0.603380 |
| Breast PCC-Spot | 0.633267 | 0.602026 |
| Breast RMSE | 0.299574 | 0.310674 |

虽然 PCC-Pro 略升，但 PCC-Spot 和 RMSE 退化，候选拒绝，不扩展七折。

运行目录：`SCP-NET-work/runs/protein_convex_simplex_anchor_Breast_Cancer_20261030/`。
