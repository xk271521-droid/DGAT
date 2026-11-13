# ProteinConvex masked reconstruction探索结果

日期：2026-10-30。用户授权尝试掩码重建以增强跨切片泛化。

## 候选A：基因空间重建

输入额外随机mask50%基因，辅助头输出32维source-fit RNA PCA分数，反变换到9764维后只在被mask基因上计算MSE，权重0.05。
20轮完成，源验证最佳PCC-Pro=0.751954（epoch20），Breast为PCC-Pro=0.584191、PCC-Spot=0.616922、RMSE=0.302345。
相比ProteinConvex基线0.602035/0.633267/0.299574，三项均退化。

该候选的目标定义存在不足：32维PCA低秩近似无法重建全部9764基因的非PCA残差，把不可重建误差混入辅助损失。因此不把它作为方法学负结果，只保留为实现审计记录。

## 候选B：标准化PCA程序重建

修正为：masked RNA输入直接预测干净source-fit PCA32程序，在标准化PCA分数空间计算MSE；PCA均值和尺度只由source-fit计算。50% gene mask、辅助权重0.05、原ProteinConvex扰动和其余训练设置不变。

20轮完成，源验证最佳为epoch19，PCC-Pro=0.752655；Breast结果为：

| 指标 | ProteinConvex基线 | masked program候选 | 候选减基线 |
|---|---:|---:|---:|
| 源验证最佳PCC-Pro | 0.751729 | 0.752655 | +0.000926 |
| Breast PCC-Pro | 0.602035 | 0.579970 | -0.022065 |
| Breast PCC-Spot | 0.633267 | 0.577747 | -0.055520 |
| Breast RMSE-Global | 0.299574 | 0.319139 | +0.019565 |

修正版仍然出现源验证略升、目标切片明显下降。说明在当前source-fit PCA32定义、50%基因遮挡、0.05辅助权重下，掩码程序重建没有改善跨切片泛化。

## 决策

两个候选均不扩展七折，不替换ProteinConvex基线。候选B是有效实现后的真实负结果；候选A只作为工程定义错误记录。
这不证明所有masked reconstruction都无效。下一次若继续，应改变辅助目标定义或加入源切片级域验证，而不是继续扫描当前权重追Breast分数。

运行目录：

- 候选A：`runs/protein_convex_masked_reconstruction_Breast_Cancer_20261030/`。
- 候选B：`runs/protein_convex_masked_program_Breast_Cancer_20261030_v2/`。
