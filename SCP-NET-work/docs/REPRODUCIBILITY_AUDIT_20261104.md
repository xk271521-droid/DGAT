# 2026-11-05 版本锁定与复现口径审计

## 当前锁定的官方版本

项目对外使用的版本名称已锁定为：

`ProteinConvex + 初始粗糙度尺度门控 + Strict RNA + Source Affine`

该名称、18 个方向的 checkpoint、Source Affine 校准产物以及当前基准指标均保留，不改动已有权重、数据集或实验参数。训练目录仍保留历史名称 `relative_roughness_gated_selection`，目录名不作为方法定义。

## 已核对的事实

- 训练权重位于 `runs/relative_roughness_gated_runs/`，共 18 个方向。
- Source Affine 校准位于 `runs/relative_roughness_gated_runs_affine_calibrated/`。
- `run_source_only_rna_affine_calibration.py` 先用源验证集拟合逐蛋白 slope/intercept，再打开目标 ADT；目标标签没有参与校准选择。
- 主训练脚本在 checkpoint 选择前只用目标 RNA/空间预测计算粗糙度，目标 ADT 在选择完成后才读取。

## 实际训练与校准逻辑

- 主训练脚本在每个 epoch 使用目标 RNA/空间预测计算无标签空间粗糙度。
- same 任务：`adaptive_alpha = 0`。
- cross 任务：若 `initial_roughness < 0.008`，则 `adaptive_alpha = 15.0`，否则为 0。
- checkpoint 选择：`selection_score = source_validation_pcc_pro - adaptive_alpha * current_absolute_roughness`。
- 目标 ADT 在 checkpoint 选择完成后才读取；Source RNA 校准的 rho、逐蛋白 lambda、slope 和 intercept 也只在源验证集上拟合。
- 目标 ADT 只用于最终指标计算，不参与轮次选择或校准参数选择。

## 当前结果与数据核验

- 当前 18 个方向的物理 checkpoint 和校准预测矩阵均已落盘。
- 独立读取预测矩阵与目标 ADT 真值后，18/18 个方向的 8 项明细指标均复核通过，最大浮点差约 `2.2e-8`。
- 18 个方向宏平均指标为：PCC-Pro `0.4694787442`、RMSE-Global `0.2885851786`、MAE-Global `0.1852496849`。
- 桌面正式表已将总览 MAE-Global 更正为 `0.1852496849`。

## 论文与复现表述

现有结果应按“初始粗糙度尺度门控”版本复现和表述。不得把这批权重和指标写成由 `Growth(e)=R(e)/R(1)` 相对增长率公式重新训练得到。若未来实现另一套相对增长率规则，应建立独立实验目录、独立结果表和独立审计记录。
