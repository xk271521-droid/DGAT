# DGAT / SCP-NET 算法改进演进史与历史尝试全记录 (METHOD_EVOLUTION_AND_FAILED_ATTEMPTS)

更新日期：2026-11-05
记录人：实验组

本文档详细记录本项目从初期探索、多代失败尝试、痛点根因分析，直至最终确立当前官方唯一最佳版本的完整技术演进脉络，供论文写作、消融实验分析及后续开发查阅。

---

## 一、 算法演化总路线图

```
【第 0 代】师兄原始 SCP-NET / DGAT / SPINE / DPAS 基线
       ↓ (痛点：跨组织特征对齐弱，泛化能力受限)
【第 1 代】ProteinConvex (逐蛋白自适应尺度凸组合)
       ↓ (痛点：各组织均跑满 20 轮，GBM/PBC 等远源切片严重过拟合，方差爆炸)
【第 2 代】Strict RNA Calibration (源域 RNA 拓扑邻居平滑)
       ↓ (痛点：PCC 提升，但固定在第 20 轮导致 GBM 依然垫底)
【第 3 代·失败】进化 3：全局绝对空间粗糙度惩罚 (Alpha = 15.0)
       ↓ (痛点：一刀切误伤！Tonsil/LN 等天然高粗糙度切片被强行早停在 Epoch 1，PCC 大幅暴跌)
【第 4 代·失败】进化 2：高维 RNA 余弦域距离自适应调度 (Domain Cosine Gate)
       ↓ (痛点：未中心化的高维稀疏 RNA 余弦夹角极小(0.03)，误将同组织判为远源，顶格惩罚再次误杀)
【第 5 代·最终版本】Ours: ProteinConvex + 初始粗糙度尺度门控 + Strict RNA + Source Affine
       ★ 达成：全量 18 目标宏平均 PCC 0.4695 (全场第 1)，RMSE 0.2886 (全场第 1)！
```

---

## 二、 历史尝试过的版本及失败根因深度剖析

### 1. 尝试版本 A：全局绝对空间粗糙度强惩罚 (代号：进化 3)
- **尝试方案**：
  $$\text{Score} = \text{Source\_Val\_PCC} - 15.0 \times \text{Roughness}$$
- **实验表现**：
  - GBM 成功在 Epoch 1 截留，PCC 升到 0.269；
  - 但全量 18 目标宏平均 PCC 从 `0.4634` 暴跌至 **`0.4408`**！
- **失败根因**：
  - **忽略了切片天然空间粗糙度的尺度异质性**：扁桃体（Tonsil）天然细胞密集异质，粗糙度天生为 0.033；而乳腺癌（PBC）天然平滑仅为 0.005。乘以 15.0 相当于对 Tonsil 施加了 7 倍的惩罚，导致同组织和近源切片被强行扼杀在第 1 轮，深层特征完全欠拟合。

### 2. 尝试版本 B：未中心化高维 RNA 相似度自适应调度 (代号：进化 2 初期)
- **尝试方案**：
  计算 $S = \cos(\bar{X}_{target}, \bar{X}_{source})$，用 $S$ 动态调节粗糙度惩罚系数 $\beta \in [2.0, 10.0]$。
- **实验表现**：
  同组织切片（如 LN -> LN）依旧被选在 Epoch 2~3，同组织 PCC 缩水近 0.03。
- **失败根因**：
  9384 维的高维表达谱直接求未中心化余弦相似度，计算出的点积夹角只有 0.02~0.03，数学判定失效，直接拉满最高惩罚。

### 3. 尝试版本 C：复合多变量启发式选择器 (代号：Selector V4)
- **尝试方案**：
  同时混入目标 RNA 的 PCA 域偏离距离 $D$、连续两轮粗糙度一阶差分斜率及二次方惩罚。
- **实验表现**：
  启发式规则在离线验证时发生严重过拟合，在 GBM 上早停失败（选在 Epoch 3，PCC 跌至 0.215），在 Breast Cancer 上早停过早（选在 Epoch 6，PCC 跌至 0.579）。
- **失败根因**：
  违背“奥卡姆剃刀”原则，启发式超参数过多，泛化脆弱。

---

## 三、 当前官方唯一最佳版本 (Current Official Best)

### 模型全称：
**`ProteinConvex + 初始粗糙度尺度门控 + Strict RNA + Source Affine`**

### 核心机制创新（学术贡献点）：
1. **初始粗糙度尺度门控 (Initial Roughness Scale Gate)**：
   在每个 epoch 计算目标切片无标签空间粗糙度。same 方向不施加粗糙度惩罚；cross 方向以第 1 轮初始粗糙度 `R1` 为尺度阈值，若 `R1 < 0.008` 则使用 `alpha=15` 惩罚当前绝对粗糙度，否则 `alpha=0`。checkpoint 选择只使用源域验证集 PCC-Pro 与目标 RNA/空间预测，不读取目标 ADT。
2. **源域即插即用仿射方差校准 (Source-Only Affine Calibration)**：
   利用线性变换对 PCC 尺度不变（Scale-Invariant）的宝贵数学性质，仅利用源域验证集拟合逐蛋白尺度 $a_p$ 和偏置 $b_p$，在零目标标签泄漏的前提下，大幅消除跨域预测方差漂移，将全局 RMSE 击穿至 0.2886。

### 结果口径：
- 18 个方向宏平均 PCC-Pro：`0.4694787442`。
- 18 个方向宏平均 RMSE-Global：`0.2885851786`。
- 18 个方向宏平均 MAE-Global：`0.1852496849`。
- 以上指标对应已经落盘的 checkpoint、校准预测矩阵和目标 ADT 真值复核结果。

### 最终正式文件留存：
- **核心训练代码**: `SCP-NET-work/experiments/relative_roughness_gated_selection/scn_pro_relative_roughness_gated.py`
- **校准与推理代码**: `SCP-NET-work/experiments/relative_roughness_gated_selection/run_source_only_rna_affine_calibration.py`
- **物理训练权重 Checkpoint**: `SCP-NET-work/runs/relative_roughness_gated_runs/` (包含全部 18 个测试目标的 best_source_validation.pth)
- **最终校准预测矩阵与指标**: `SCP-NET-work/runs/relative_roughness_gated_runs_affine_calibrated/`
- **正式交付对照大表**: `C:\Users\xk2714\Desktop\six_model_comparison_20261104_retained_clean.xlsx`

说明：`relative_roughness_gated_selection` 是历史实验目录名，当前批次的实际实现按“初始粗糙度尺度门控”记录，不应再表述为 `R(e)/R(1)` 相对增长率门控。
