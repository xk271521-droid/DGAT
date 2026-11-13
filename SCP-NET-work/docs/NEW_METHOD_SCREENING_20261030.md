# 新方法筛选记录（2026-10-30）

## 当前判断

ProteinConvex 的主要瓶颈仍是跨切片泛化：源验证集可以升高，但完整未见目标切片不一定同步升高。已有记录已经排除了延长 epoch、RNA 程序重建、MetaLODO、source-conditional PCC、MaskedCognate、CenterAdapter、MixStyle/Fishr、均衡损失/采样和若干输出/隐藏一致性路线。

gene-mask reconstruction 候选的第一版得到 Breast `PCC-Pro=0.584191`，但实现把32维 PCA 重建先反变换到全基因后逐基因计算损失，混入了低秩近似误差；第二版又没有标准化 PCA 分数，辅助损失量级压制了蛋白主损失。因此这条路线目前不能作为干净的负结果，本记录不把它当作成功或失败。

## 优先候选 A：NicheTrans 风格的 RNA 噪声强度对照

官方 NicheTrans 的 Breast 参数使用 `noise_rate=0.2`；其输入层通过 element-wise dropout 扰动 RNA。当前 ProteinConvex 使用 `gene_dropout_rate=0.5`，且每次训练前向都启用。此前只做过同时关闭 RNA dropout 和邻居 mask 的联合对照，不能推出中间噪声强度无效。

第一轮只改一个参数：

```text
gene_dropout_rate: 0.5 -> 0.2
neighbor_mask_rate: 0.25（保持）
内部 dropout: 0.3（保持）
```

数据、seed=42、20 epoch、训练/验证拆分、source-only checkpoint 选择和目标测试流程全部保持不变。先跑 Breast_Cancer 单折，三项指标均不退化才考虑七折。

实际结果：源验证最佳 `0.744234`，Breast `PCC-Pro/PCC-Spot/RMSE=0.569231/0.564962/0.316891`，三项均低于 ProteinConvex，已拒绝，不扩展七折。详见 `RESULT_NICHE_NOISE02_BREAST_20261030.md`。

静态尺度 simplex 约束也已完成：Breast `0.603380/0.602026/0.310674`，PCC-Spot 和 RMSE 退化，已拒绝。详见 `RESULT_SIMPLEX_ANCHOR_BREAST_20261030.md`。

DPAS 风格 SmoothL1 也已完成：Breast `0.590217/0.606240/0.302931`，三项退化，已拒绝。源验证宏平均选点最终仍选第20轮，目标结果与基线完全一致，未产生变化。

DPAS 风格 RNA feature-KNN regional 邻居也已完成：Breast `0.599391/0.588623/0.309864`，三项退化，已拒绝。

## 优先候选 B：NicheTrans 风格的邻居腐蚀

当前 ProteinConvex 对每个 batch 的邻居 token 独立做 `p=0.25` attention mask，但 CognateAnchor 的 local/regional RNA 均值仍使用未遮挡邻居。NicheTrans 的 Breast 代码是约30%的 batch 触发邻居腐蚀，触发时每个邻居向量以 Bernoulli(0.5) 整体置零。

若候选 A 无效，再做候选 B：同一邻居 mask 同时用于 attention 和 CognateAnchor 均值，保持 `gene_dropout_rate=0.5`。这不是完整复现 NicheTrans，只是来源清楚的单变量扰动策略对照。

## 中优先候选 C：DPAS-Graph 双图 RNA 相似邻居

当前 ProteinConvex 只使用空间 KNN。DPAS-Graph 还构造 RNA PCA 特征 KNN（PCA 保留约85%方差、k=10，并使用边权），再让模型学习空间关系和表达关系的贡献。可在 source-fit RNA PCA 上构图，目标只使用 RNA，不读取目标 ADT；通过零初始化的有界 gate 接入现有局部/区域分支。

该方案改动较大，涉及图构造和显存，排在候选 A/B 后面。

## 证据位置

- NicheTrans Breast 参数：`upstream/NicheTrans/args/args_breast_cancer.py`
- NicheTrans RNA 噪声：`upstream/NicheTrans/model/nicheTrans_hd.py`
- NicheTrans 邻居腐蚀：`upstream/NicheTrans/utils/utils_training_breast_cancer.py`
- 当前 ProteinConvex 扰动入口：`SCP-NET-work/experiments/protein_convex/scn_pro_protein_convex.py`
- 当前额外扰动关闭后的负结果：`SCP-NET-work/docs/NO_EXTRA_MASK_BREAST_20261030.md`

本记录只提出候选和受控实验协议，不把任何候选描述为已经提升分数。
