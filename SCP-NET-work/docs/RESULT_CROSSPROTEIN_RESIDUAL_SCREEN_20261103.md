# 低秩蛋白间残差解码器结构筛选：负结果

日期：2026-11-03

## 候选结构

在每个蛋白独立预测头后增加 rank-8 低秩蛋白间残差：

```text
prediction + W_up(GELU(W_down(prediction)))
```

`W_up` 零初始化，因此初始化时输出与原始 ProteinConvex 结构相同。候选参数增加 496 个。原正式脚本、模型、权重和结果均未覆盖。

## 固定筛选协议

- 面板：9384 RNA 基因、31 ADT 蛋白；
- seed 42、20 epochs、batch 64、validation ratio 0.1、lr 0.0005、aux weight 0.25、PCC weight 0.05，其余设置与正式 ProteinConvex 一致；
- source-only validation 选 checkpoint；target ADT 只用于最终预测评分；
- 困难方向：`cross_holdout_Breast_Cancer -> PBC_PR_6837` 和 `cross_holdout_GBM -> GBM`；
- 两方向都在训练后直接用验证过的校准程序进行 source-only 校准前，本记录只评估了结构候选的未校准输出。

## 结果

| 目标 | 指标 | 未校准基线 | 残差结构候选 | 变化 |
|---|---|---:|---:|---:|
| PBC_PR_6837 | PCC-Pro | 0.243771 | 0.227997 | -0.015774 |
| PBC_PR_6837 | RMSE-Global | 0.289749 | 0.287940 | -0.001810 |
| PBC_PR_6837 | MAE-Global | 0.180434 | 0.179873 | -0.000560 |
| GBM | PCC-Pro | 0.203952 | 0.188468 | -0.015484 |
| GBM | RMSE-Global | 0.317649 | 0.319601 | +0.001952 |
| GBM | MAE-Global | 0.209802 | 0.211427 | +0.001625 |

PBC 上只有误差指标轻微改善，主要相关性退化；GBM 的 PCC-Pro、RMSE 和 MAE 均退化。候选未达到两个困难方向都改善的筛选条件，判为**未通过，不扩展全 18 目标**。

说明：此表与当前 source-only RNA 校准基线不是完全同一推理 pipeline；候选没有再做推理校准，因为结构本身已在两个目标的 PCC-Pro 上均低于未校准基线。不能把它报告为校准组合候选的结论。

## 清理

按用户要求提交了删除候选脚本和运行产物的操作，但自动审核拒绝了对以下候选专属目录的递归删除，因此截至本记录更新时间，代码、日志、checkpoint、预测和配置仍保留在原位：

- `SCP-NET-work/experiments/protein_convex_crossprotein_residual_20261103/`
- `SCP-NET-work/runs/protein_convex_crossprotein_residual_20261103/`

负结果记录和项目索引保留；正式 baseline、校准结果及五模型/校准版工作簿未改动。候选产物不能描述为已删除。
