# ProteinConvex + Strict RNA Calibration：跨组织测试集最佳 epoch 诊断

本实验只针对当前最佳模型 `ProteinConvex + Strict RNA Calibration`，完整覆盖跨组织 10 个目标：LN 2 个、Tonsil 4 个、Breast Cancer 3 个、GBM 1 个。

每个方向仍按正式 ProteinConvex 配置训练 20 个 epoch。每个 epoch 都重新用 source validation 自动选择 RNA/空间校准参数，再用目标 ADT 计算该 epoch 的校准后 RMSE、MAE 和 PCC。目标 ADT 因此参与了 epoch 选择，结果是 test-set oracle 诊断，不是正式公平结果。

正式基线 checkpoint 仍按 source validation PCC-Pro 选择；本目录只回答：如果为了 RMSE/MAE 事后选择测试集最佳 epoch，当前模型是否有明显潜在改善。
