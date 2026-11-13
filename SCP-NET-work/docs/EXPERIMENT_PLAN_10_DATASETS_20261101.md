# 十个数据集同组织与跨组织实验规划

日期：2026-11-01

## 当前目标

使用 10 个数据集，比较当前最佳方法 `ProteinConvex + Strict RNA Calibration` 与四个作者方法：

- DGAT
- DPAS-Graph
- NicheTrans
- SPINE

四个作者方法用于公平比较分数，不替代当前主方法。所有实验必须在同一数据注册表、同一测试分组、同一蛋白面板和统一指标下记录。

## 统一数据面板

10 个数据集为：

- `A1_LN`
- `D1_LN`
- `A1_TNSL`
- `D1_TNSL`
- `HumanTonsil_10x`
- `HumanTonsil_AddOns_S1`
- `Breast_Cancer`
- `PBC_PR_6837`
- `PBC_PR_6835_5A`
- `GBM`

10 个数据集共同 RNA 面板为 **9384 个基因**，共同 ADT 面板为 **31 个蛋白**。原来 7 个主数据集的 9764 基因面板不再用于这套 10 数据集正式实验。

由于输入基因维度发生变化，不能直接加载旧的 9764 基因 checkpoint。五个模型都必须按照各自官方实现和当前固定实验协议重新训练；这属于重新训练/重新适配输入维度，不是把旧 checkpoint 继续微调后冒充新结果。

## 同组织实验

按实验规划文档规定的四组，每组做双向训练和测试：

| 组别 | 方向一 | 方向二 |
|---|---|---|
| LN | `A1_LN` 训练 → `D1_LN` 测试 | `D1_LN` 训练 → `A1_LN` 测试 |
| TNSL | `A1_TNSL` 训练 → `D1_TNSL` 测试 | `D1_TNSL` 训练 → `A1_TNSL` 测试 |
| Human Tonsil | `HumanTonsil_10x` 训练 → `HumanTonsil_AddOns_S1` 测试 | `HumanTonsil_AddOns_S1` 训练 → `HumanTonsil_10x` 测试 |
| Breast Cancer | `PBC_PR_6837` 训练 → `PBC_PR_6835_5A` 测试 | `PBC_PR_6835_5A` 训练 → `PBC_PR_6837` 测试 |

`Breast_Cancer` 和 `GBM` 在当前 10 个数据集中各只有一个样本，因此不单独组成同组织配对。

同组织实验共 8 个方向；5 个模型全部运行，共 **40 个模型训练方向**。

## 跨组织实验

按组织类别分为四组，每次完整留出一个组织类别：

| 留出组织 | 测试集 | 训练集 |
|---|---|---|
| LN | `A1_LN`、`D1_LN` | 其余 8 个数据集 |
| Tonsil | `A1_TNSL`、`D1_TNSL`、`HumanTonsil_10x`、`HumanTonsil_AddOns_S1` | 其余 6 个数据集 |
| Breast Cancer | `Breast_Cancer`、`PBC_PR_6837`、`PBC_PR_6835_5A` | 其余 7 个数据集 |
| GBM | `GBM` | 其余 9 个数据集 |

每个留出方向只用训练组织拟合预处理统计量、选择 checkpoint 和选择模型参数。测试组织的 ADT 只能用于最终指标计算。

跨组织实验共 4 个留出方向；5 个模型全部运行，共 **20 个模型训练方向**。

同组织和跨组织合计 **60 个模型训练方向**。每个方向都必须保存数据分组、输入面板、参数、checkpoint、日志和指标。

## 当前本地可运行性盘点

本机硬件：NVIDIA GeForce RTX 5060 Laptop GPU，约 8 GB 显存；当前 `st2p` 环境为 Python 3.9.23、PyTorch 2.8.0+cu128。

| 模型 | 源码状态 | 当前环境状态 | 本地结论 |
|---|---|---|---|
| ProteinConvex + Strict RNA | 项目内已有训练脚本 | 当前 `st2p` 环境可用 | 可作为本地优先模型，但要先切换到 9384 基因输入并重新训练 |
| DGAT | 官方源码已归档 | 需要 `torch_geometric` 等包，当前环境未安装；作者代码还检查多 GPU 分支 | 不能直接启动，需独立环境和输入适配 |
| DPAS-Graph | 官方源码已归档 | 官方环境是 Linux/Python 3.9/旧 CUDA 组合，当前 Windows 环境缺少 `torch_geometric`、`torch_scatter`、`torch_sparse` | 不在当前环境直接跑，优先服务器独立环境 |
| NicheTrans | 官方源码已归档 | 依赖较多，作者参数文件带有作者机器路径，输入格式不是当前统一 h5ad 直接格式 | 可做独立适配后再评估，本地暂不启动 |
| SPINE | 官方源码已归档 | README 推荐 Python 3.12.9、PyTorch 2.7；默认 500 epoch、batch size 1，要求 RNA embedding 数据入口 | 本地只适合后续小规模 smoke test；正式 10 数据集长任务放服务器 |

## 当前明确不做的事情

- 不直接使用旧的 9764 基因 checkpoint；
- 不把外部数据补零到 9764 维后混入新结果；
- 不修改四个作者仓库的原始源码；
- 不在当前环境里强行安装一套混杂依赖后直接跑正式结果；
- 不在未完成输入适配和 smoke test 前启动 60 个正式训练方向。

## 下一步顺序

1. 建立 10 数据集统一注册表、9384 基因顺序和 31 蛋白顺序；
2. 为五个模型分别建立独立适配目录，不改官方源码归档；
3. 先对每个模型做单方向小规模 smoke test，检查输入形状、显存、输出和指标；
4. 确认本地可运行模型与服务器模型的分工；
5. 再启动同组织 8 方向和跨组织 4 方向的正式训练。

2026-11-01：已启动本地 `ProteinConvex + Strict RNA` 十数据集版本。输入为 9384 个共同 RNA 基因和 31 个共同 ADT 蛋白，运行目录为 `SCP-NET-work/runs/protein_convex_10_dataset_20261101/`，日志目录为 `SCP-NET-work/experiments/protein_convex_10_dataset_20261101/logs/`。当前第一方向为 `A1_LN` 训练到 `D1_LN` 测试，后台运行正常。

2026-11-01：NicheTrans 官方模型已完成 1 epoch smoke test，确认 9384→31 输入、空间 12 邻居和作者模型前向可运行；smoke 产物已删除。正式 12 个方向已在 `SCP-NET-work/runs/nichetrans_10_dataset_20261101/` 后台启动。使用作者乳腺癌参数：noise_rate=0.2、dropout_rate=0.2、40 epoch、batch size=32、Adam、lr=0.0003、weight_decay=0.0005、seed=1。官方源码未修改，适配代码在 `SCP-NET-work/experiments/nichetrans_10_dataset_20261101/`。

2026-11-01：DPAS-Graph 适配已完成 1 epoch smoke test。官方双图模型、图构建和损失保持不变；只适配 9384/31 共同面板、十数据集方向、源数据 10% validation 选点和结果保存。正式运行目录为 `SCP-NET-work/runs/dpas_graph_10_dataset_20261101/`，使用官方默认 100 epochs、`lr_scale=1.0`、`seed=0`；目标 ADT 延迟到 checkpoint 选定后读取。适配代码在 `SCP-NET-work/experiments/dpas_graph_10_dataset_20261101/`，官方源码未修改。

2026-11-01：本地仅完成 `same_LN_A1_to_D1` 后暂停。剩余 11 个 DPAS-Graph 方向已转到服务器 `/public/liqiang/DGAT` 排队运行，使用单卡物理 GPU 0；GPU 2、3 不使用，方向按顺序串行执行。服务器结果目录为 `/public/liqiang/DGAT/SCP-NET-work/runs/dpas_graph_10_dataset_20261101_server/`。

2026-11-02：服务器检查确认 DPAS-Graph 剩余 11 个方向全部完成并保存结果；所有方向日志均有 `COMPLETE`，没有发现错误。连同本地完成的 `same_LN_A1_to_D1`，DPAS-Graph 12 个方向全部完成。

2026-11-03：SPINE 的 `cross_holdout_Breast_Cancer` 初版训练虽已完成，但旧编排把两个 PBC 数据集误列为 source，未按乳腺癌组织整体留出。已按正确协议在服务器 v6 重新训练：7 个数据集训练，`Breast_Cancer`、`PBC_PR_6837`、`PBC_PR_6835_5A` 三个数据集测试。训练和评估均完成，三组预测与真值已拉回本地 `results/SPINE_10_dataset_20261101/analysis_ready/cross_holdout_Breast_Cancer/`；旧 v2/v3/v5 中间目录和错误结果已清理。

2026-11-03：在保持 ProteinConvex checkpoint、模型结构、训练参数、数据和 9384/31 面板不变的情况下，完成十数据集 source-only RNA 邻居推理校准。每个方向只用 source validation 选择图混合参数及逐蛋白平滑强度；目标 ADT 只在参数选定后最终评估。基线重建与保存的 18 个正式预测最大绝对差异为 9.54e-7。校准后 18/18 目标的 PCC-Pro、RMSE-Global、MAE-Global 均改善；完整记录见 `SCP-NET-work/experiments/protein_convex_rna_calibration_10_dataset_20261103/RESULTS.md`。校准预测独立保存于 `SCP-NET-work/runs/protein_convex_rna_calibration_10_dataset_20261103/all_directions/`，未校准基线保持不变。推荐比较工作簿为 `results/five_models_20261103/comparison/five_model_comparison_20261103_rna_calibrated.xlsx`，原始五模型比较表保留作基线。
