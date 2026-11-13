# 蛋白等权方差标准化损失：Breast_Cancer

本候选基于 ProteinConvex，仅把最终监督项改为：

```text
0.5 * 原始 MSE + 0.5 * 按源训练集蛋白标准差标准化后的 MSE
```

标准差只由 6 个源切片的 fit spots 计算。模型结构、辅助 MSE、PCC 项、seed=42、20 epoch、mask 和 source-only checkpoint 选择均保持不变。

| 指标 | ProteinConvex 基线 | 蛋白等权标准化损失 |
|---|---:|---:|
| 源验证最佳 PCC-Pro | 0.751729 | 0.764238 |
| Breast PCC-Pro | 0.602035 | 0.589329 |
| Breast PCC-Spot | 0.633267 | 0.601366 |
| Breast RMSE | 0.299574 | 0.304175 |

源验证分数提高，但目标切片三项指标均退化，候选拒绝，不扩展七折。该结果说明当前训练目标与 PCC-Pro 存在错位，但简单的 50/50 方差标准化没有解决跨切片泛化问题。

运行目录已按用户要求清理；原始 ProteinConvex 基准未修改。
