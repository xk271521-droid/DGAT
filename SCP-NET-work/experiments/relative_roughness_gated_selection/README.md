# 初始粗糙度尺度门控自适应选择模型 (Initial Roughness Scale-Gated Selection)

日期：2026-11-05

## 一、 背景与前序实验问题
在“进化3”及“宏观域相似度”方案中，我们发现采用绝对粗糙度惩罚存在严重的组织异质性偏差：
- 扁桃体（Tonsil）和淋巴结（LN）等结构丰富的近源组织天生具备较高的空间粗糙度（0.015 ~ 0.033），而乳腺癌（PBC）初始粗糙度仅为 0.005；
- 采用绝对数值惩罚公式 $\text{Score} = \text{PCC} - \beta \times \text{Roughness}$ 时，无论 $\beta=15.0$ 还是 $\beta=10.0$，都会将天生高粗糙度的近源组织一刀切强行扼杀在第 1~2 轮，导致同质组织 PCC 整体大幅暴跌 0.03~0.05。

## 二、 实际执行公式与无作弊自适应机制
代码路径：`experiments/relative_roughness_gated_selection/scn_pro_relative_roughness_gated.py`

1. **第 1 轮初始尺度**：
   在第 1 轮记录目标切片无标签空间粗糙度 `initial_roughness = Roughness(1)`。

2. **初始粗糙度尺度门控**：
   - **same 任务**：`adaptive_alpha = 0`，只使用源域验证集 PCC-Pro 选择 checkpoint；
   - **cross 任务**：当 `initial_roughness < 0.008` 时使用 `adaptive_alpha = 15.0`，否则使用 `adaptive_alpha = 0`；
   - 当前选择分数为：
   $$\text{Penalty}(e) = \alpha \times \text{Roughness}(e)$$
   $$\text{Selection Score}(e) = \text{Source\_Val\_PCC}(e) - \text{Penalty}(e)$$

   该规则是当前 18 个方向 checkpoint 和指标的实际生成逻辑。目录名中的 `relative_roughness_gated` 是历史名称，不代表本批产物使用了相对增长率公式。

3. **零数据泄漏保证**：
   整个训练与轮次选择过程中，目标切片的 ADT 蛋白标签绝对处于封闭状态，仅在选定最终 checkpoint 之后进行一次性推理评估。Strict RNA 和 Source Affine 参数也只在源域验证集上拟合。
