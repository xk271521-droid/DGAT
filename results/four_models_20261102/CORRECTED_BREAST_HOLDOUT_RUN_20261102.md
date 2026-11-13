# Breast Cancer 严格留出修正运行记录

日期：2026-11-02

## 固定协议

- source：`A1_LN`、`D1_LN`、`A1_TNSL`、`D1_TNSL`、`HumanTonsil_10x`、`HumanTonsil_AddOns_S1`、`GBM`
- targets：`Breast_Cancer`、`PBC_PR_6837`、`PBC_PR_6835_5A`
- 不把任何一个 PBC/乳腺癌数据放进 source；同一个 checkpoint 评估该组三个 targets。

## NicheTrans

- 状态：完成
- 运行目录：`SCP-NET-work/runs/nichetrans_10_dataset_20261102_corrected_breast_holdout/cross_holdout_Breast_Cancer`
- source/target：与上述固定协议一致
- 输入面板：9384 RNA 基因、31 ADT 蛋白
- 参数：作者配置保持不变，40 epochs、batch size 32、neighbors 12、noise/dropout 0.2、lr 0.0003、weight decay 0.0005、seed 1
- best epoch：35
- source-validation PCC-Pro：0.7293701
- 三个目标 prediction/truth 已保存，形状分别为 4160×31、4850×31、3243×31

## DPAS-Graph

- 状态：服务器运行中
- 结果目录：`/public/liqiang/DGAT/SCP-NET-work/runs/dpas_graph_10_dataset_20261102_corrected_breast_holdout`
- 脚本：`experiments/dpas_graph_10_dataset_20261101/run_corrected_breast_holdout_server_20261102.sh`
- source/target：服务器适配器 `cross_holdout_Breast_Cancer` 方向已修正并确认 source 为上述 7 个数据集
- 参数：官方 DPAS-Graph 模型/损失保留，100 epochs、lr_scale 1.0、seed 0，单卡 GPU 0
- 尚未完成前不把它的旧泄漏结果或新部分日志纳入正式表格

2026-11-02 启动状态复核：修正版适配器已上传并确认 source 列表为 7 个正确数据集。为减少重复预处理，复用了哈希匹配的六个 source fit/validation 图缓存；GBM 的新 split 图按修正版 seed 重新构建。最后一次 SSH 状态检查时为 89/100 epochs，source-validation PCC-Pro=0.6350，checkpoint 已保存，仍需完成训练、3 个 target 预测/指标与结果落盘。按最近训练速度估计剩余约 10–30 分钟（含最终 target 图构建和预测），该估计有运行时波动。SSH 已断开，服务器后台任务继续运行。

## 旧结果处理

旧目录中的两个错误结果已按用户要求删除。本地 NicheTrans 旧目录 `SCP-NET-work/runs/nichetrans_10_dataset_20261101_fixed_panel/cross_holdout_Breast_Cancer/`、本地 DPAS-Graph 旧目录 `SCP-NET-work/runs/dpas_graph_10_dataset_20261101_all_results/cross_holdout_Breast_Cancer/` 以及服务器 DPAS-Graph 旧目录 `/public/liqiang/DGAT/SCP-NET-work/runs/dpas_graph_10_dataset_20261101_server/cross_holdout_Breast_Cancer/` 已移除。旧指标没有保留在四模型表，相关 6 行已清空并标记 `已删除 / 待重跑`。
