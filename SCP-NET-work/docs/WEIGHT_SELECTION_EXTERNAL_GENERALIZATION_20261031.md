# 外部泛化测试权重选择记录（2026-10-31）

> 历史记录：本权重选择对应旧的 9764 基因映射方案，不适用于新的 10 数据集共同面板实验。

适用方法：`SCN-Pro-ProteinConvex + Strict RNA Calibration`

## 一、神经网络 checkpoint

外部测试统一使用 `Breast_Cancer` 留一折训练得到的冻结 checkpoint：

`SCP-NET-work/runs/protein_convex_Breast_Cancer/best_source_validation.pth`

该 checkpoint 的训练源切片为：

- `A1_LN`
- `D1_LN`
- `A1_TNSL`
- `D1_TNSL`
- `HumanTonsil_10x`
- `HumanTonsil_AddOns_S1`

checkpoint 只按这 6 个源切片的 validation spots 平均 PCC-Pro 选择：

- 最佳 epoch：`20`
- 源验证平均 PCC-Pro：`0.751729`

`PBC_PR_6837`、`PBC_PR_6835_5A` 和 `GBM` 没有参与 checkpoint 选择或重新训练。

## 二、Strict RNA 校准权重

校准参数只使用上述 6 个源切片的 validation spots 选择。

### RNA/空间图混合权重

搜索网格：

```text
rho ∈ {0, 0.25, 0.5, 0.75, 1.0}
```

最终选择：

```text
rho = 0
```

含义：最终使用纯 RNA 表达相似图，不混入空间距离图权重。

### 逐蛋白扩散强度

每个蛋白独立搜索：

```text
lambda_p ∈ {0, 0.1, 0.2, 0.3, 0.4, 0.5}
```

最终选择的逐蛋白 `lambda_p` 保存在：

`SCP-NET-work/experiments/external_generalization_strict_calibration_20261031/selected_lambda_per_protein.npy`

平均值：

```text
mean(lambda_p) = 0.396774
```

### 校准公式

```text
calibrated_prediction_p
  = (1 - lambda_p) × raw_prediction_p
  + lambda_p × RNA_similarity_neighbour_prediction_p
```

## 三、外部数据使用边界

三个外部数据集只使用 RNA 和空间坐标进行前向预测、切片上下文计算和 RNA 图构建。外部 ADT 在预测和校准完成后才读取，仅用于最终指标计算。

外部数据不参与：

- 模型训练；
- checkpoint/epoch 选择；
- `rho` 选择；
- `lambda_p` 选择；
- 任何超参数调优。

## 四、复现文件

- 评估脚本：`SCP-NET-work/experiments/external_generalization_strict_calibration_20261031/evaluate_external_strict_calibration.py`
- 选择网格：`SCP-NET-work/experiments/external_generalization_strict_calibration_20261031/source_selection_grid.csv`
- 逐蛋白权重：`SCP-NET-work/experiments/external_generalization_strict_calibration_20261031/selected_lambda_per_protein.npy`
- 外部结果：`SCP-NET-work/docs/RESULT_EXTERNAL_GENERALIZATION_STRICT_CALIBRATION_20261031.md`
