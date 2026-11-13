# 旧外部映射清理记录

日期：2026-11-01

## 清理内容

已删除旧外部推理方案中由“映射到 9764 基因面板并补零”产生的派生内容：

- `experiments/external_generalization_strict_calibration_20261031/`
- `experiments/pbc_external_validation_20261030/external_inference.py`
- `experiments/pbc_external_validation_20261030/results/`

## 保留内容

- 原始 PBC 和 GBM 数据；
- PBC 和 GBM 的独立预处理矩阵；
- 原始数据构建和预处理脚本；
- 旧结果文档作为历史审计记录，但不再作为新的 10 数据集正式结果。

## 新输入协议

10 个数据集统一使用 9384 个共同 RNA 基因。新的模型输入、训练、验证、测试和汇报结果必须基于这一共同面板重新生成。
