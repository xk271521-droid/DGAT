# 无效候选清理记录

日期：2026-10-30

根据用户要求，删除已确认没有提升效果的候选实验代码和运行产物，保留结果记录文档、官方源码、原始 ProteinConvex 代码和正式基准结果。

已删除的候选包括：NicheTrans 噪声 `gene_dropout=0.2`、静态尺度 simplex、DPAS SmoothL1、源切片宏平均选点、DPAS RNA feature-KNN regional 邻居，以及 masked reconstruction / masked program 候选。

保留的记录文档包括 `RESULT_*.md`、`NEW_METHOD_SCREENING_20261030.md` 和本文件。原始基准目录 `SCP-NET-work/runs/protein_convex_Breast_Cancer/` 未修改。

## 追加清理（2026-10-30）

应用户要求继续删除已确认无提升的候选代码和运行产物：

- `experiments/protein_convex_target_masked_20261030/` 及对应运行目录（目标掩码适应，Breast PCC-Pro −0.0433）。
- `experiments/protein_convex_balanced_20261030/` 及对应运行目录（逐蛋白均衡 MSE，PCC-Pro −0.0172）。
- `runs/protein_convex_Breast_Cancer_50epoch_20261030/`（延长训练，PCC-Pro −0.0137）。
- `runs/protein_convex_no_extra_mask_Breast_Cancer_20261030/`（关闭额外扰动，PCC-Pro −0.0711）。
- 上述候选日志，以及此前残留的 masked program / masked reconstruction 日志。

保留：失败记录 `RESULT_*.md`、正向的 `experiments/spatial_smoothing_20261030/`、五模型原始代码与运行目录。`experiments/protein_convex_dann_20261030/` 尚未评估，不属于“已确认无效”，本次未删除。
