# DGAT 十数据集服务器实验

本目录使用 `upstream/DGAT` 中固定提交的官方 `GATEncoder`、`Decoder_mRNA`、`Decoder_Protein` 和图构建代码；官方源码保持原样。

必要适配只有：固定 9384 基因/31 蛋白面板、当前十数据集方向、单 GPU 运行、源数据 validation 选 checkpoint，以及保存最终蛋白预测和真值。训练损失权重、隐藏层、dropout、优化器和 100 epoch 默认值沿用作者代码。

## 2026-11-02 结果审计

服务器目录 `runs/dgat_10_dataset_20261101_server/` 的 12 个方向均已完成，保存了 18 个目标数据集的预测和真值文件，且没有残留 DGAT 训练进程。

审计发现 `cross_holdout_Breast_Cancer` 的方向配置把 `PBC_PR_6837` 和 `PBC_PR_6835_5A` 列入了 sources，同时又把它们列入 targets。按照组织级留出协议，乳腺癌组织的三个数据集必须全部留出，正确 sources 应为其余 7 个数据集。因此该方向虽已训练完成，但存在乳腺癌组织内数据泄漏，暂不纳入正式公平比较；其余 11 个方向的 source/target 划分符合当前实验计划。

已启动独立修正版重跑，输出目录为 `runs/dgat_10_dataset_20261101_server_corrected_breast/`。修正版 sources 为 `A1_LN`、`D1_LN`、`A1_TNSL`、`D1_TNSL`、`HumanTonsil_10x`、`HumanTonsil_AddOns_S1`、`GBM`，targets 为 `Breast_Cancer`、`PBC_PR_6837`、`PBC_PR_6835_5A`；仍使用官方 100 epochs、seed=0、9384/31 面板和单 GPU 运行。

该次首次启动因 GPU0 被其他任务占满，在第一个 batch 发生 CUDA out-of-memory，未生成结果；错误输出保留在 `runs/dgat_10_dataset_20261101_server_corrected_breast/launcher.log`。随后已改用空闲 GPU2，在 `runs/dgat_10_dataset_20261101_server_corrected_breast_retry_gpu2/` 重新启动，source/target 配置和模型参数不变。
