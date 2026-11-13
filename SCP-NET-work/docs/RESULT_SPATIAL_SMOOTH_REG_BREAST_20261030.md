# 训练期空间预测平滑正则：Breast_Cancer

本候选基于 ProteinConvex，在源 fit spots 的真实空间邻居边上增加预测一致性 MSE，权重固定为 `0.05`。目标 RNA/ADT 不参与该正则项，其他结构、损失、seed、20 epoch、mask 和 source-only checkpoint 选择保持不变。

| 指标 | ProteinConvex 基线 | 训练期空间平滑正则 |
|---|---:|---:|
| 源验证最佳 PCC-Pro | 0.751729 | 0.754244 |
| Breast PCC-Pro | 0.602035 | 0.602284 |
| Breast PCC-Spot | 0.633267 | 0.595128 |
| Breast RMSE | 0.299574 | 0.313460 |

PCC-Pro 仅增加 `0.000249`，但 PCC-Spot 和 RMSE 明显退化，候选拒绝，不扩展七折。此前推理期空间平滑的正向结果没有直接转化为训练期一致性正则的收益。

候选代码和运行目录已清理，原始 ProteinConvex 基准未修改。
