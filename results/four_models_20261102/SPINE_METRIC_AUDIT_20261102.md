# SPINE 同组织 PBC 指标复核

审计日期：2026-11-02

对 `same_PBC_6837_to_6835` 和 `same_PBC_6835_to_6837`，分别核对了 SPINE 的 `*_metrics.json`、`summary.json` 以及保存的 prediction/truth 数组。

| 方向 | PCC-Pro | PCC-Spot | Spearman-Pro | RMSE-Global |
|---|---:|---:|---:|---:|
| `PBC_PR_6837 -> PBC_PR_6835_5A` | -0.0093 | -0.0369 | -0.0042 | 0.5377 |
| `PBC_PR_6835_5A -> PBC_PR_6837` | 0.1532 | 0.1281 | 0.1400 | 0.2588 |

负值是有效的 Pearson/Spearman 相关结果。第一条方向的 31 个蛋白平均相关接近零且略偏负，PCC-Spot 也略偏负，说明预测与目标之间没有稳定线性关系；这不是指标符号写反或表格错位。该方向的误差指标也明显更高，说明确实是一次失败的跨切片泛化，而不是单个相关系数偶然为负。

SPINE 的 checkpoint 按 source validation 选择，目标 ADT 只用于最终评估；因此该结果不能通过目标标签选点解释，也不应被删除或改成零。
