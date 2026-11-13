# DPAS-Graph 十数据集实验

本目录是 `upstream/DPAS-Graph` 官方实现的独立适配层。官方源码保持原样；适配只负责读取项目统一输入、建立实验方向、源验证选点和结果保存。

## 固定协议

- 10 个数据集共同 RNA 面板：9384 个基因；
- 共同 ADT 面板：31 个蛋白；
- 同组织 8 个双向方向，跨组织 4 个留出方向；
- 每个方向从源数据中固定抽取 10% validation spots；
- checkpoint 只按源 validation 的 mean PCC-Pro 选择，目标 ADT 只在最终预测后读取；
- 官方 DPAS-Graph 模型、双图构建、损失权重和优化器参数保持不变；
- 每个方向保存 checkpoint、训练历史、最终蛋白预测和真值、逐蛋白/逐 spot 指标。

运行环境为项目内 `.venv-dpas`，基于现有 `st2p` 的 PyTorch 2.8.0+cu128，并补充与其匹配的官方 PyG wheel。

2026-11-01：`same_LN_A1_to_D1` 1 epoch smoke test 已通过。官方双图构建、GPU 前向、源验证 PCC-Pro 选点和目标预测文件保存均正常。正式运行保持官方默认 100 epochs、`lr_scale=1.0`、`seed=0`；目标 ADT 在源验证 checkpoint 选定后才读取。首次正式启动在第 5 epoch 被人工停止，未作为结果使用；后续使用新输出目录重新启动。

2026-11-01：正式本地运行已按要求暂停。有效完成结果仅为 `same_LN_A1_to_D1`，最佳 epoch=94，源验证 PCC-Pro=0.494359，输入为 9384 基因/31 蛋白。`same_LN_D1_to_A1` 目录只包含被停止后的不完整运行，不纳入结果。服务器准备包已上传到 `/public/liqiang/DGAT/server_runs/dpas_graph_20261101`，服务器任务尚未启动。

服务器排队协议：只使用一张 GPU（当前指定 GPU 0，不使用 GPU 2、3），剩余方向按脚本顺序串行执行；允许与其他任务共享该 GPU，但本 DPAS-Graph 进程本身不申请多卡。

2026-11-01：服务器队列已启动。已完成的 `same_LN_A1_to_D1` 不再重复运行；剩余 11 个方向从 `same_LN_D1_to_A1` 开始，输出目录为 `/public/liqiang/DGAT/SCP-NET-work/runs/dpas_graph_10_dataset_20261101_server/`。服务器使用 `HyperST` 环境和 `CUDA_VISIBLE_DEVICES=0`，当前无 GPU 2/3 占用。

2026-11-02：服务器队列检查确认剩余 11 个方向全部生成 `results.json`，11 个 `run.log` 均有 `COMPLETE` 标记，DPAS-Graph 进程已退出；未发现 traceback/error。加上本地已完成的 `same_LN_A1_to_D1`，DPAS-Graph 十数据集协议的 12 个方向已全部完成。

2026-11-02 protocol audit：原 `cross_holdout_Breast_Cancer` 把两个 PBC 数据集错误地同时列入 sources 与 targets，三个目标的旧结果不用于公平比较。严格修正版正在服务器独立目录 `/public/liqiang/DGAT/SCP-NET-work/runs/dpas_graph_10_dataset_20261102_corrected_breast_holdout/` 运行，direction 为 7 个正确 source（LN、Tonsil、GBM）到 3 个 Breast Cancer targets。官方模型、损失和参数保持不变，使用 100 epochs、lr_scale=1.0、seed=0、单卡 GPU 0。服务器适配器脚本为 `run_corrected_breast_holdout_server_20261102.sh`，适配方向列表也已修正，避免后续误启动旧 source 划分。
