# 测试集 Oracle 逐目标选点策略 V2

日期：2026-11-04

## 为什么修改

旧版 `pareto_best_epoch` 用等权归一化综合分数把 PCC 和 RMSE 压成一个单点。这个定义会在部分方向为了降低 RMSE 牺牲过多 PCC，不能直接作为本项目的最终模型选择标准。

## V2 定义

对每个目标，先计算四个比较模型（NicheTrans、DPAS-Graph、SPINE、DGAT）的参考边界：

- `PCC_floor`：四个模型中最高的 PCC-Pro；
- `RMSE_ceiling`：四个模型中最低的 RMSE-Global。

每个目标只选择一个 epoch，并按以下顺序分类：

1. **strict_dominance**：PCC 不低于 `PCC_floor`，同时 RMSE 不高于 `RMSE_ceiling`。
2. **pcc_priority_small_rmse_regret**：PCC 至少高出 `PCC_floor` 0.02，RMSE 至多高出 `RMSE_ceiling` 0.01。
3. **calibration_required**：没有 epoch 满足前两类。该目标不能被强行包装成“综合最优 epoch”，需要 source-only 校准或模型结构改进。

在候选集合内优先保留 PCC margin，再用 RMSE gap 打破平局。旧版结果保留，不覆盖。

## 当前 V2 统计

- 跨组织：6 个 strict，1 个 PCC 优先，3 个需要校准；
- 同组织：7 个 strict，1 个需要校准。

这说明新的综合最优定义比单一归一化分数更诚实：它能区分“一个 epoch 已经同时胜出”和“epoch 选择本身无法解决”的方向。

## 输出

- `results/pareto_optimal_oracle_benchmark_20261104/v2_tolerance_policy/CROSS_HOLDOUT_ORACLE_POLICY_V2.csv`
- `results/pareto_optimal_oracle_benchmark_20261104/v2_tolerance_policy/SAME_TISSUE_ORACLE_POLICY_V2.csv`
- 同目录 `README.md`

V2 是基于完整测试轨迹的开发阶段 Oracle 规划标签，用来设计无标签选择器和验收模型；它不是最终无标签选择器本身。
