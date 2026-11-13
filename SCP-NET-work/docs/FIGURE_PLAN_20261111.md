# 最终模型绘图计划与口径记录

更新日期：2026-11-11

本文档只记录当前正式结果的绘图方案，不修改训练代码、模型参数、数据集或既有结果。

## 当前正式方法

正式主方法固定为：

`ProteinConvex + 初始粗糙度尺度门控 + Strict RNA + Source Affine`

核心训练代码：

`experiments/relative_roughness_gated_selection/scn_pro_relative_roughness_gated.py`

Source Affine 校准代码：

`experiments/relative_roughness_gated_selection/run_source_only_rna_affine_calibration.py`

最终训练产物位于：

`runs/relative_roughness_gated_runs/`

最终校准预测和指标位于：

`runs/relative_roughness_gated_runs_affine_calibrated/`

正式评估口径为 10 个数据集、9384 个共同 RNA 基因、31 个共同 ADT 蛋白，共 18 个目标方向：8 个同组织方向和 10 个跨组织方向。目标 ADT 只能在 checkpoint、邻域校准参数和逐蛋白仿射参数全部由 source validation 冻结后读取。

## 代码链路摘要

1. 每个 source section 按 seed 打乱并划分 fit/validation；target RNA 可用于无标签推理，target ADT 延迟到最后评估。
2. 每个 spot 构造 8 个无自环邻居：前 4 个为 local，后 4 个为 regional，并给邻居加入按中位最近邻距离归一化的相对位置。
3. center、local、regional 三路 RNA 经过共享 RNA encoder；local/regional 经过 protein-conditioned multi-head attention。
4. global router 加低秩 protein-specific residual 产生每个蛋白的三尺度 softmax 权重，融合 center/local/regional latent。
5. RNA-only slice context 来自 source-fit RNA PCA 的均值、标准差和三个 RNA/geometry 尺度描述；它进入 context encoder、router 和低秩 context adapter。
6. bounded slice canonicalization 对每片 RNA 做受限残差校正；cognate RNA anchor 对可匹配的蛋白注入 latent 先验及直接预测。
7. 训练损失为 final MSE + 0.25×三路 auxiliary MSE + 0.05×protein-wise PCC loss；训练时使用 feature dropout、neighbor token mask 和梯度裁剪。
8. checkpoint 选择分数为 source-validation PCC-Pro 减去无标签 target spatial roughness penalty。same 方向 alpha=0；cross 方向只有初始 roughness < 0.008 时 alpha=15，否则 alpha=0。这里是“初始粗糙度尺度门控”，不是相对增长率 `R(e)/R(1)`。
9. Source Affine 在 source validation 上搜索空间/RNA 邻域混合系数 rho、逐蛋白平滑系数 lambda，并拟合逐蛋白 slope/intercept；冻结后才应用于 target。

## 建议的主图

### Figure 1：方法结构图

画一张 forward DAG，而不是把所有类名堆在一起：

`center RNA + local RNA + regional RNA + relative position`
→ `shared RNA encoder`
→ `protein-conditioned local/regional attention`
→ `center/local/regional level states`
→ `protein-specific 3-scale router`
→ `cognate RNA anchor + protein-specific MLP`
→ `31-protein prediction`。

从侧边接入 `RNA-only slice context`，标明它同时作用于 router 和 context adapter；训练外层再画 `roughness-gated source-validation checkpoint selection` 和 `source-only affine calibration`。该图解释结构，不能把后处理画成网络内部可学习层。

### Figure 2：18 个目标的正式基准

主面板使用最终校准结果，按 same/cross 分面：

- 每个目标用点或短横线展示 PCC-Pro、RMSE-Global、MAE-Global；
- 旁边给出六模型宏平均点图，模型包括最终版本、消融基线、NicheTrans、DPAS-Graph、SPINE、DGAT；
- 目标顺序固定为实验协议顺序，不能按结果排序后隐藏任务结构。

主结果的已核对宏平均为：全部 18 个目标 PCC-Pro=0.4694787442、RMSE-Global=0.2885851786、MAE-Global=0.1852496849；同组织 8 个目标 PCC-Pro=0.5340850725、RMSE-Global=0.2598233763；跨组织 10 个目标 PCC-Pro=0.4177936817、RMSE-Global=0.3115946203。

### Figure 3：最终版本相对消融的增益

用 paired slope/dumbbell 或 delta heatmap 对比 `relative_roughness_gated_runs` 的未校准预测和 `relative_roughness_gated_runs_affine_calibrated` 的最终预测。三张小面板分别展示 ΔPCC-Pro、ΔRMSE-Global、ΔMAE-Global；误差指标的改善方向要在图例中明确为“越低越好”。

### Figure 4：初始粗糙度尺度门控

第一面板画 18 个方向的 `initial_roughness`，在 0.008 处画阈值线，并用颜色标示 same/cross 和 alpha 是否开启。真实产物中触发 alpha=15 的是三个跨组织方向：PBC_PR_6835_5A、PBC_PR_6837、GBM，三者 best epoch 都为 1；same 方向即使 roughness 较低也不触发门控。

第二面板选 GBM、PBC 和一个 LN/TNSL 方向，画 epoch 对应的 source-validation PCC、target roughness、selection score，并标出 source-selected best epoch。不要画成 target ADT 参与选择的曲线。

### Figure 5：Source Affine 校准

第一面板展示每个方向 baseline→calibrated 的 PCC/RMSE/MAE paired change；第二面板展示 selected rho；第三面板展示 31 个蛋白的 lambda heatmap 或分布。图注写明所有参数来自 source validation，target ADT 只用于最终评分。

### Figure 6：三尺度路由与空间定性结果

- `protein_scale_attention` 可画成 target direction × protein 的 center/local/regional 热图，或按组织类别画三尺度均值分布；
- local/regional attention 的原始形状为 `[spots, heads, proteins, neighbors]`，先在 spots/head 上聚合，再画 protein × neighbor；对接近精确均匀的 attention 计算熵或标准差并过滤，不能默认解释为生物学重要性；
- 选 GBM、Breast_Cancer、HumanTonsil_AddOns_S1 各一个目标，读取 RNA h5ad 的真实 `obsm['spatial']`，画 truth、prediction、residual 的空间图，再配 2–4 个代表蛋白的 truth-vs-prediction 散点。

## 建议放入补充材料的图

- 31 个蛋白的 PCC/RMSE/MAE 全量热图；
- 18 个方向的完整 metrics_history 小 multiples；
- 每个方向的 best epoch、source validation 与 target-oracle 轨迹对照；
- `rna_slice_context_descriptors.npz` 和 canonicalization 的诊断图；
- 全部模型在 same/cross/all 三种分组下的 8 项指标表。

## 不能直接画或不能直接解释的内容

1. `results/pareto_optimal_oracle_benchmark_20261104/` 中的 target-oracle 轨迹只能标为开发期诊断或上限，不能当成正式无标签选择结果。
2. `results.json` 中的 `wall_seconds` 是 `time.monotonic()` 的绝对读数，不能直接当训练耗时。
3. 目标空间图必须使用真实坐标；不能把 `[n_spots, 31]` 预测矩阵伪装成空间坐标。
4. 当前最终 attention 文件没有保存 `slice_context`、`anchor_gate` 或 center/local/regional 分支 prediction。若后续需要这些量，必须写独立只读 extraction 脚本，用冻结 checkpoint 重新 forward，并把派生产物单独记录。
5. 当前校准脚本正文的默认 `RUN_ROOT/OUTPUT_ROOT` 仍指向旧的 ProteinConvex calibration 目录，而最终产物位于 `relative_roughness_gated_runs*`。正式画图和复现前必须以最终目录中的 `calibration_results.json`、配置和脚本快照核对路径，不能直接无检查地运行默认命令。

## 推荐执行顺序

先做 Figure 2、Figure 3 和 Figure 4，确认主结果、增益和门控逻辑的口径；再做 Figure 5 的校准解释；最后做 Figure 1 和 Figure 6 的结构/定性图。所有图都从现有结果派生，不重新训练，不改变参数，不打开新的 target 标签。
