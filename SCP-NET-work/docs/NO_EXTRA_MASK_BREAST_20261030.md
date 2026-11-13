# ProteinConvex关闭额外RNA/邻居扰动对照

日期：2026-10-30。用户授权验证当前项目额外训练扰动是否影响Breast泛化。

## 实验定义

从ProteinConvex原始Breast_Cancer 20 epoch运行独立派生，唯一改变：

- `gene_dropout_rate=0.5 -> 0`
- `neighbor_mask_rate=0.25 -> 0`

保留网络内部`dropout=0.3`。网络结构、ProteinConvex alpha、数据、7切片注册表、
6源1目标、源fit/validation拆分、seed42、batch64、学习率、损失权重、20 epoch和
source-validation PCC-Pro checkpoint选点全部不变。目标ADT只在选点后读取。

这不是DGAT或DPAS-Graph官方代码复现，只是判断当前项目额外扰动是否影响ProteinConvex跨切片结果的受控对照。

## 预期判读

关闭扰动后如果目标测试提升，说明当前额外增强可能与本项目输入/任务不匹配；
如果不提升，问题更可能在跨切片分布和source-validation代理指标，而不是这两项mask。
该单折结果不作为正式七折结论，只有确认有价值后才考虑完整七折。

运行目录：`SCP-NET-work/runs/protein_convex_no_extra_mask_Breast_Cancer_20261030/`。

## 实际结果

| 指标 | 原ProteinConvex（0.5/0.25） | 关闭额外扰动 | 候选减基线 |
|---|---:|---:|---:|
| 源验证最佳PCC-Pro | 0.751729（epoch20） | 0.741112（epoch13） | -0.010617 |
| Breast PCC-Pro | 0.602035 | 0.530923 | -0.071111 |
| Breast PCC-Spot | 0.633267 | 0.561866 | -0.071402 |
| Breast RMSE-Global | 0.299574 | 0.321843 | +0.022269 |
| Breast MAE-Global | 0.202944 | 0.218243 | +0.015299 |

20轮正常完成，候选由源验证集选择第13轮。基线和候选的seed、20 epoch、
模型参数量、数据拆分数组和评估预测形状均核对一致；候选唯一配置差异是两项额外扰动均为0。

## 决策

关闭额外RNA feature dropout和邻居token mask后，源验证和Breast测试均明显变差。
因此在当前ProteinConvex和当前数据协议下，这两项项目级训练增强不能删掉。
该结果不证明它们是DGAT或DPAS-Graph作者的原始设置，也不证明50%一定是最优值；
它只说明“直接关闭两项”不是有效改进。该候选不扩展七折，不替换ProteinConvex基线。

运行产物已按 2026-10-30 清理要求删除；本文保留结果记录。
