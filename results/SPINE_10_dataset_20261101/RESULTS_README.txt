SPINE 10数据集结果（本地副本）

来源：服务器 /public/liqiang/DGAT/SCP-NET-work/runs/spine_10_dataset_20261101_v2,v3,v5，以及修正版 v6
本地目录：E:\DGAT\results\SPINE_10_dataset_20261101
分析目录包含 12 个实验方向；乳腺癌跨组织方向已由 v6 修正版补齐 3 个目标数据集。

说明：
- `analysis_ready/` 是后续分析的唯一入口，包含同组织 8 个方向和跨组织 4 个方向的最终预测、真值和指标。
- `cross_holdout_Breast_Cancer` 已按组织级留出规则包含 `Breast_Cancer`、`PBC_PR_6837` 和 `PBC_PR_6835_5A` 三个目标。
- `audit/cross_holdout_Breast_Cancer_v6/` 保存修正版的配置、数据清单和训练/评估日志。
