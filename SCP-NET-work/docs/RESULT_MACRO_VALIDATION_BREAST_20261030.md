# 源切片宏平均验证选点：Breast_Cancer

模型和训练完全保持 ProteinConvex 基线，只把 checkpoint 选择分数改为 6 个源切片验证集 PCC-Pro 的宏平均，避免按 spot 混池时大切片权重更高。

| 指标 | ProteinConvex 基线 | 源切片宏平均选点 |
|---|---:|---:|
| 选点源验证分数 | 0.751729（混池） | 0.721695（宏平均） |
| Breast PCC-Pro | 0.602035 | 0.602035 |
| Breast PCC-Spot | 0.633267 | 0.633267 |
| Breast RMSE | 0.299574 | 0.299574 |

宏平均方案最终仍选择第20轮，因此目标预测与基线完全一致。该方案没有带来变化，不扩展七折。

运行目录：`SCP-NET-work/runs/protein_convex_macro_validation_Breast_Cancer_20261030/`。
