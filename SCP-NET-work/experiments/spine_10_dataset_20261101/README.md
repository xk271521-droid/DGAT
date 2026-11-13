# SPINE 十数据集适配

官方 SPINE 源码位于 `upstream/SPINE`，本目录只保存数据适配和实验编排代码，不修改官方模型源码。

适配将统一的 9384 RNA 基因和 31 ADT 蛋白 h5ad 转换为官方需要的 `RNA_EMBED/*.h5`、`protein.h5ad`、基因/蛋白列表和 `train_0.csv` / `test_0.csv`。

官方训练入口使用 `--normalize_method None`，避免对已冻结的 DGAT 预处理矩阵重复 log1p。

每个 held-out 样本的预测保存在模型目录下的 `heldout_predictions/`，包括 `.npy` 矩阵、带 spot barcode 和蛋白名的 `.h5ad`、压缩 CSV，以及真值和指标文件。训练队列结束后，`save_missing_predictions_server.sh` 会在 CPU 上串行补齐任何缺失的预测，避免占用训练 GPU。

## 2026-11-02 组织级留出复核

服务器 `spine_10_dataset_20261101_v5` 的 `cross_holdout_Breast_Cancer` 训练已经完成，但旧编排脚本误把 `PBC_PR_6837` 和 `PBC_PR_6835_5A` 同时列入 source。它们因此被写成 `__fit/__val`，没有生成真正的 held-out 目标输入；该方向存在乳腺癌组织内数据泄漏，不纳入正式比较。

已保留 v5 目录作为审计记录，并按正确协议启动独立的 `spine_10_dataset_20261101_v6`：乳腺癌组织的 3 个数据集全部留出，剩余 7 个数据集训练，目标为 `Breast_Cancer`、`PBC_PR_6837` 和 `PBC_PR_6835_5A`。模型参数、9384 基因/31 蛋白面板和训练流程保持不变。
