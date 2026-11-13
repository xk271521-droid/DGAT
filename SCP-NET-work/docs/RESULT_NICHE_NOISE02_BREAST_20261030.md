# NicheTrans 风格 RNA 噪声强度对照：Breast_Cancer

日期：2026-10-30

## 唯一改动

在 ProteinConvex 原始 Breast_Cancer 外层留一折配置中，仅将：

```text
gene_dropout_rate: 0.5 -> 0.2
```

保持 `neighbor_mask_rate=0.25`、内部 `dropout=0.3`、seed=42、20 epoch、数据拆分、损失、source-only checkpoint 选择和冻结预处理不变。

## 结果

| 指标 | ProteinConvex 基线 | gene dropout=0.2 |
|---|---:|---:|
| 源验证最佳 PCC-Pro | 0.751729 | 0.744234 |
| Breast PCC-Pro | 0.602035 | 0.569231 |
| Breast PCC-Spot | 0.633267 | 0.564962 |
| Breast RMSE | 0.299574 | 0.316891 |

候选在第19轮达到源验证最佳，但目标切片三项指标均退化，因此拒绝该设置，不扩展七折。

运行目录：`SCP-NET-work/runs/protein_convex_niche_noise02_Breast_Cancer_20261030/`。
