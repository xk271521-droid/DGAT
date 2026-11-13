# 跨组织测试集最佳 epoch 诊断：已完成部分

日期：2026-11-03

## 范围与状态

本诊断只针对 `ProteinConvex + Strict RNA Calibration`。截至中断，完整跑完 3/10 个跨组织目标：`A1_LN`、`D1_LN`、`A1_TNSL`。`D1_TNSL` 方向仅有初始化产物，没有 epoch 日志或 `results.json`，不计入结果。其余 6 个跨组织目标未运行。

每个已完成方向均按原设置训练 20 个 epoch；每个 epoch 都重新由 source validation 选择 RNA/空间校准参数，再在 target 上记录校准后的指标。测试集 ADT 被用于逐 epoch 评分和事后选择，因此以下是 oracle 诊断，不是正式泛化结果。

## RMSE / MAE 对照

作者模型列为正式六模型工作簿的 NicheTrans、DPAS-Graph、SPINE、DGAT。表中的“作者最佳”是这四者中相应误差最低者。

| 目标 | source validation 选中 epoch | 校准后 RMSE / MAE | 测试集最佳 RMSE（epoch） | 测试集最佳 MAE（epoch） | 四作者最佳 RMSE | 四作者最佳 MAE |
|---|---:|---:|---:|---:|---:|---:|
| A1_LN | 18 | 0.340914 / 0.214195 | 0.337667（2） | 0.214195（18） | 0.347890（NicheTrans） | 0.218869（NicheTrans） |
| D1_LN | 18 | 0.295455 / 0.187224 | 0.295455（18） | 0.187224（18） | 0.291111（NicheTrans） | 0.188372（NicheTrans） |
| A1_TNSL | 20 | 0.353886 / 0.227629 | 0.336772（1） | 0.214959（2） | 0.342080（NicheTrans） | 0.213923（NicheTrans） |

前三个目标等权均值：source-validation 校准结果 RMSE `0.330085`、MAE `0.209683`；逐目标 test-oracle 最佳 RMSE `0.323298`、最佳 MAE `0.205459`。这是测试标签选择后的诊断上限，不能与正式模型结果等同。

## 初步解读

- `A1_LN` 和 `A1_TNSL` 的测试集最低 RMSE 分别在 epoch 2 和 epoch 1，早于 source-validation PCC 选出的 epoch 18 和 20；误差分别改善约 `0.00325` 和 `0.01711`。
- `D1_LN` 的 RMSE/MAE 最佳 epoch 正好是 source-selected epoch 18，没有因测试集选点获益。
- 与四个作者模型比较，测试集 oracle RMSE 在 3 个目标中赢过作者最佳者 2 个；MAE 在 3 个中赢过 2 个。`D1_LN` 的 RMSE 仍不及 NicheTrans，`A1_TNSL` 的 MAE 仍略高于 NicheTrans。
- 结果初步支持“按 PCC 选择 epoch 可能错过较低 RMSE 的 checkpoint”，但目前只有 3/10 跨组织目标，不能下总体结论。

## 可复核文件

逐 epoch 指标、source-selected checkpoint 和 oracle 预测分别保存在 `SCP-NET-work/runs/protein_convex_calibrated_testset_oracle_20261103_cross/` 各方向目录。对照作者模型的来源是 `results/six_models_20261103/comparison/six_model_comparison_20261103.xlsx`。
