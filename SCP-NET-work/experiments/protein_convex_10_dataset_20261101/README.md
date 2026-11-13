# ProteinConvex 十数据集运行记录

本目录运行 `ProteinConvex + Strict RNA` 的 10 数据集版本。

输入固定为 9384 个十数据集共同 RNA 基因和 31 个共同 ADT 蛋白。

训练脚本每个方向会保存：

- `best_source_validation.pth`：源验证集选择的 checkpoint；
- `results.json`：训练、选择和审计信息；
- `target_source_selected_predictions.npy`：源验证选定 checkpoint 对测试代表切片产生的最终原始蛋白预测；
- `target_source_selected_per_protein_metrics.csv`：逐蛋白指标；
- `target_source_selected_attention.npz`：注意力和尺度信息；
- `script_snapshot.py`、`config.json`、`source_split_indices.npz`：复现记录。

这些预测文件必须保留。后续 Strict RNA Calibration 和统一结果分析直接读取已保存预测，不重复训练或重新推理。

正式方向包括 8 个同组织双向方向和 4 个跨组织留出方向，共 12 个 ProteinConvex 训练方向。
