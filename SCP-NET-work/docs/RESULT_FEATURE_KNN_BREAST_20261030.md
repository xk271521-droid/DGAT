# DPAS 风格 RNA 相似邻居：Breast_Cancer

保留每个 spot 的4个空间局部邻居，把4个区域邻居改为同一 RNA panel 上的表达相似邻居。目标切片只用 RNA 建图；模型、损失、seed、20 epoch 和 source-only checkpoint 选择不变。

| 指标 | ProteinConvex 基线 | RNA feature-KNN regional |
|---|---:|---:|
| 源验证最佳 PCC-Pro | 0.751729 | 0.746714 |
| Breast PCC-Pro | 0.602035 | 0.599391 |
| Breast PCC-Spot | 0.633267 | 0.588623 |
| Breast RMSE | 0.299574 | 0.309864 |

三项目标指标均退化，候选拒绝，不扩展七折。

运行目录：`SCP-NET-work/runs/protein_convex_feature_knn_Breast_Cancer_20261030/`。
