# ProteinConvex 多专家 RNA/空间推理融合：七折探索

本候选在原始 ProteinConvex 预测、空间邻居平滑预测和 RNA 相似邻居预测之间，为每个蛋白分别选择两个融合权重。参数只用当前外层实验的六个源切片 validation spots 选择。

| 指标 | 多专家融合 | 严格 RNA 校准 |
|---|---:|---:|
| 七折 PCC-Pro | 0.594753 | 0.594800 |
| 七折 PCC-Spot | 0.644613 | 0.644555 |
| 七折 RMSE-Global | 0.271663 | 0.272659 |
| Breast PCC-Pro | 0.612526 | 0.612301 |
| Breast PCC-Spot | 0.638632 | 0.642322 |
| Breast RMSE-Global | 0.290239 | 0.290963 |

多专家融合只在 RMSE 和 Breast PCC-Pro 上略有优势，没有同时超过严格 RNA 校准的相关性指标，因此当前主结果继续使用严格 RNA 校准。多专家代码已清理，逐折 CSV 保留在本实验目录中供审计。
