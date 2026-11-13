# ProteinConvex + 源切片域对抗：Breast_Cancer

本候选基于 `SCN-Pro-ProteinConvex`，只新增源切片域分类器和梯度反转层，域损失权重 `0.05`、GRL 系数 `1.0`。原 ProteinConvex 的结构、扰动、损失、seed=42、20 epoch、source-only checkpoint 选择和冻结数据均保持不变。

| 指标 | ProteinConvex 基线 | ProteinConvex + DANN |
|---|---:|---:|
| 源验证最佳 PCC-Pro | 0.751729 | 0.748694 |
| Breast PCC-Pro | 0.602035 | 0.589101 |
| Breast PCC-Spot | 0.633267 | 0.586269 |
| Breast RMSE | 0.299574 | 0.320748 |

三项目标指标均退化，拒绝该候选，不扩展七折。该结果说明当前源切片标签的全局域对抗约束可能抹掉了对蛋白预测有用的真实组织差异。

运行目录：`SCP-NET-work/runs/protein_convex_dann_Breast_Cancer_20261030/`。
