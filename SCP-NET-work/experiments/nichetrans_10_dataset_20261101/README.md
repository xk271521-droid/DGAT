# NicheTrans 十数据集实验

官方模型代码从 `upstream/NicheTrans/model/nicheTrans.py` 原样导入。

保留的作者参数：`noise_rate=0.2`、`dropout_rate=0.2`、40 epoch、batch size 32、Adam、learning rate 0.0003、weight decay 0.0005、seed 1、空间 KNN 12 邻居。

必要适配：输入维度改为固定共同面板中的 9384 个 RNA 基因和 31 个 ADT 蛋白；从当前统一 h5ad 读取；按空间坐标构建 12 邻居；使用源切片 10% validation 的 PCC-Pro 选择 checkpoint；保存所有目标切片的预测和真值。基因和蛋白顺序直接读取 `SCP-NET-work/docs/ten_dataset_common_panel_20261101.json`。

官方源码未修改。运行脚本只负责数据适配和实验编排。

2026-11-01：首次重跑的前几个方向启动早于固定面板代码更新，已停止且不作为正式结果。当前正式重跑目录为 `SCP-NET-work/runs/nichetrans_10_dataset_20261101_fixed_panel/`，12 个方向全部从头使用同一份 9384 基因/31 蛋白固定面板。

2026-11-02 protocol audit：原 `cross_holdout_Breast_Cancer` 把 `PBC_PR_6837` 和 `PBC_PR_6835_5A` 错误地同时放入 sources 与 targets，结果不用于公平比较。严格修正版已使用独立目录 `SCP-NET-work/runs/nichetrans_10_dataset_20261102_corrected_breast_holdout/cross_holdout_Breast_Cancer/` 完成，sources 为 LN 2 个、Tonsil 4 个和 GBM 共 7 个数据集；targets 为 `Breast_Cancer`、`PBC_PR_6837`、`PBC_PR_6835_5A`。作者模型参数未改，fixed 9384/31 面板、40 epochs、seed 1。该目录为正式 Breast Cancer leave-one-tissue 结果。
