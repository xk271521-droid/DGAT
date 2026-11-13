# DGAT / SCP-NET 最终项目计划与核心记录规范 (PROJECT_PLAN)

更新日期：2026-11-05
基准配置：10 个空间转录组数据集，共同 RNA 面板 9384 基因，共同 ADT 31 蛋白。

---

## 一、 当前项目官方唯一主方法 (Official Champion Model)

本项目的主方法已正式固定并升级为：
**`ProteinConvex + 初始粗糙度尺度门控 + Strict RNA + Source Affine`**

- **模型地位**：作为项目对外汇报、论文主实验（Main Benchmark）展示的最终唯一优胜模型。
- **基线消融模型**：原版 `ProteinConvex + Strict RNA Calibration`（即未加相对粗糙度门控和未加源域仿射校准的原始版本）严格保留在对比大表中，作为消融实验（Ablation Study）对照组，证明自适应门控与尺度对齐的独立增益。
- **四个公开对比模型**：严格保留 `NicheTrans`、`DPAS-Graph`、`SPINE`、`DGAT` 官方完整对比结果。

---

## 二、 核心实验指标与验收口径

全量 18 个测试目标切片（8 个同组织 + 10 个跨组织）的最终官方表现如下：

1. **宏观综合对比**：
   - 18 目标宏平均 PCC-Pro: **0.4695**（全场最高，彻底突破旧版瓶颈）
   - 18 目标全局 RMSE-Global: **0.2886**（全场最低残差）
   - 18 目标宏平均 MAE-Global: **0.1852**
2. **跨组织留出（10 个目标切片）**：
   - 宏平均 PCC-Pro: **0.4178**（大幅碾压对手近 0.05）
   - 全局 RMSE-Global: **0.3116**（反超对手）
   - 攻克了 GBM（PCC 0.2688 / RMSE 0.2720）和 PBC_6837（PCC 0.2901 / RMSE 0.2451）的早期过拟合问题。
3. **同组织近源（8 个目标切片）**：
   - 宏平均 PCC-Pro: **0.5341**，全局 RMSE-Global: **0.2598**（8/8 全胜第一）。

---

## 三、 文件索引与维护规范

1. **桌面正式大表**:
   `C:\Users\xk2714\Desktop\six_model_comparison_20261104_retained_clean.xlsx`
2. **物理产物位置**:
   - 权重: `SCP-NET-work/runs/relative_roughness_gated_runs/`
   - 预测矩阵与校准指标: `SCP-NET-work/runs/relative_roughness_gated_runs_affine_calibrated/`
   - 核心代码: `SCP-NET-work/experiments/relative_roughness_gated_selection/`
3. **技术演进与失败方案记录**:
   详见 `SCP-NET-work/docs/METHOD_EVOLUTION_AND_FAILED_ATTEMPTS.md`。

4. **复现口径审计**:
   `SCP-NET-work/docs/REPRODUCIBILITY_AUDIT_20261104.md` 记录了当前实际门控逻辑、目标标签隔离和指标口径。现有权重和指标按实际初始粗糙度尺度门控版本使用，不得表述为相对增长率公式重新训练得到。
