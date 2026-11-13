# DPAS 风格 SmoothL1 蛋白损失：Breast_Cancer

只将最终蛋白预测项从 MSE 改为 `SmoothL1Loss(beta=1.0)`；辅助 MSE、PCC 项、模型结构、数据、seed=42、20 epoch 和 source-only 选点保持不变。

| 指标 | ProteinConvex 基线 | SmoothL1 |
|---|---:|---:|
| 源验证最佳 PCC-Pro | 0.751729 | 0.749533 |
| Breast PCC-Pro | 0.602035 | 0.590217 |
| Breast PCC-Spot | 0.633267 | 0.606240 |
| Breast RMSE | 0.299574 | 0.302931 |

三项目标指标均退化，候选拒绝，不扩展七折。

运行目录：`SCP-NET-work/runs/protein_convex_huber_Breast_Cancer_20261030/`。
