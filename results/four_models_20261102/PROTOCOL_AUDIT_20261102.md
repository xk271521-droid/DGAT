# 四模型实验协议审计

审计日期：2026-11-02

## 发现的问题

NicheTrans 和 DPAS-Graph 的 `cross_holdout_Breast_Cancer` 结果配置错误：

- targets：`Breast_Cancer`、`PBC_PR_6837`、`PBC_PR_6835_5A`
- 实际 sources：包含 `PBC_PR_6837`、`PBC_PR_6835_5A`
- 正确 sources：
  `A1_LN`、`D1_LN`、`A1_TNSL`、`D1_TNSL`、`HumanTonsil_10x`、`HumanTonsil_AddOns_S1`、`GBM`

这造成了 Breast Cancer 组织内部数据泄漏，涉及四模型表中的 6 条明细：两个模型 × 三个目标。当前 NicheTrans 和 DPAS-Graph 的这 6 条结果标记为无效，不进入最终公平比较。

## 已核对为正确的部分

- ProteinConvex 的 Breast Cancer 留出三个方向均使用 7 个正确 source；
- SPINE 的 v6 Breast Cancer manifest 明确使用同样的 7 个 source 和 3 个 targets；
- NicheTrans 和 DPAS-Graph 的 LN、Tonsil、GBM source/target 列表符合正式计划；
- 四模型表的明细数值可以从本地 prediction/truth 复算，但这不能弥补错误 source 配置造成的公平性问题。

## 当前表格状态（2026-11-03）

NicheTrans 的三个乳腺癌目标已由严格 7-source/3-target 修正版替换。DGAT 结果已拉回并通过 source/target、数组形状和指标复算审计。

当前比较表纳入 ProteinConvex + Strict RNA、NicheTrans、SPINE、DGAT；每个模型 18 个目标，合计 72 条有效明细。DPAS-Graph 按用户要求暂不纳入本版。旧错误 NicheTrans/DPAS-Graph Breast Cancer 结果已删除，不再留在表内。
