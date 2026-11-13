# 十数据集 source-only RNA 邻居校准实验计划

日期：2026-11-03

## 目标

在冻结的十数据集 ProteinConvex/SCN-Pro checkpoint 上，测试 source-only 选择的 RNA 相似邻居校准是否改善跨切片/跨组织表达预测。已有七数据集严格 source-validation 校准曾得到正结果；本计划是在当前 9384/31 新面板和正式 12 方向上重新验证，不预设一定提升。

## 固定项

- 数据：当前十数据集共同处理版本，不增删或重预处理数据；
- 面板：9384 RNA 基因、31 ADT 蛋白，沿用每个正式方向的顺序；
- 模型结构、checkpoint、训练过程、训练参数、seed、source split 全部冻结；
- target 只使用 RNA 和空间坐标建立前向输入/邻居图；target ADT 在校准参数选定后才读取，仅用于最终评分；
- 每个方向单独用其 source validation 选择参数，禁止跨方向读取 target 指标调参。

## 候选方法

每个 slice 从 RNA PCA 距离和空间距离分别构造 RBF 邻居权重。source-validation 网格选择图混合比例 `rho ∈ {0, 0.25, 0.5, 0.75, 1}`，并为每个蛋白选择 `lambda_p ∈ {0, 0.1, 0.2, 0.3, 0.4, 0.5}`。预测为原始 ProteinConvex 输出与邻居加权输出的凸组合。网格继承项目既有七数据集严格校准实现，不在本轮扩展搜索。

## 验证顺序

1. 用已保存 checkpoint、配置和数据重新生成一个方向的原始 target prediction；与正式运行保存的 `target_source_selected_predictions.npy` 比较，最大绝对差异须小于 `2e-5`。
2. 复现检查通过后，对 18 个目标逐方向校准；保存 source-validation 选参、基线/校准预测和最终指标。
3. 分别汇总 8 个同组织目标、10 个跨组织目标，以及四个跨组织留出组。主要指标为 PCC-Pro、PCC-Spot、Spearman-Pro、Spearman-Spot；RMSE/MAE 为误差护栏。
4. 先作为独立候选报告；只有相对冻结基线在预先定义的相关性汇总上提高且 RMSE 不明显退化，才讨论是否更新主方法。

## 风险与限制

邻居平滑可能模糊边界、损害 spot 级指标；source validation 的空间/表达邻居关系不保证完全代表 unseen target。结果只用于当前十数据集协议，不外推为普遍改进结论。
