# ProteinConvex 测试集最佳 epoch 诊断

日期：2026-11-03

本实验应师兄建议建立，用来回答一个诊断问题：如果每个 epoch 都直接查看目标切片 ADT，并事后选择目标集 RMSE/MAE/PCC 最好的 epoch，ProteinConvex 的误差指标理论上能达到什么水平。

这不是正式公平评估。目标 ADT 在训练过程中被读取并参与每个 epoch 的评分，因此结果属于 test-set oracle / upper-bound diagnostic，不能替换正式表格，也不能作为论文主结果。

实验完整覆盖 12 个方向和 18 个目标：8 个同组织双向方向，以及 LN、Tonsil、Breast Cancer、GBM 四个整体留出方向。训练数据、9384 基因/31 蛋白面板、seed=42、20 epochs、batch size=64、source split、优化器和模型参数均沿用正式 ProteinConvex 配置。

正式基线仍使用 source validation 的 mean PCC-Pro 选择 checkpoint。每个方向的结果会分别记录 source-validation checkpoint、目标集 PCC 最佳 epoch、RMSE 最佳 epoch 和 MAE 最佳 epoch，并保存对应预测矩阵。
