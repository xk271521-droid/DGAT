# 外部泛化测试：ProteinConvex + Strict RNA Calibration

> 历史记录：本结果使用旧的 9764 基因外部映射方案。2026-11-01 起不再作为 10 数据集正式结果。

日期：2026-10-31

## 协议

- 使用 `Breast_Cancer` 留一折训练得到的冻结 ProteinConvex checkpoint；
- 校准参数只由该 checkpoint 的六个源切片 validation spots 选择；
- 外部数据只使用 RNA 和空间坐标生成前向预测、RNA 表达图和校准结果；
- 外部 ADT 在预测和校准完成后才读取，仅用于最终指标计算；
- `PBC_PR_6837`、`PBC_PR_6835_5A` 和 `GBM` 不参与训练、checkpoint 选择或参数调优。

源验证选择结果：`rho=0`，平均 `lambda_p=0.396774`。

完整权重选择规则见：`SCP-NET-work/docs/WEIGHT_SELECTION_EXTERNAL_GENERALIZATION_20261031.md`。

## 结果

| 数据集 | spots | 公共基因面板覆盖 | 原始 PCC-Pro | 校准后 PCC-Pro | 原始 PCC-Spot | 校准后 PCC-Spot | 原始 RMSE | 校准后 RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `PBC_PR_6837` | 4,850 | 9,477/9,764 | 0.225725 | **0.236671** | 0.189028 | **0.206154** | 0.289713 | **0.267830** |
| `PBC_PR_6835_5A` | 3,243 | 9,595/9,764 | 0.246840 | **0.257423** | 0.253184 | **0.268870** | 0.368745 | **0.356226** |
| `GBM` | 5,637 | 9,594/9,764 | 0.126247 | **0.129087** | 0.076215 | **0.078898** | 0.354192 | **0.341094** |

三个外部数据集宏平均：

```text
校准后 PCC-Pro       0.207727
校准后 RMSE-Global   0.321717
```

## 解释

1. Strict RNA Calibration 在三个外部数据集上都提高了 PCC-Pro，并降低了全局 RMSE。
2. 两个新乳腺癌病人的结果高于 GBM，说明跨病人泛化明显好于跨组织泛化。
3. GBM 的绝对相关性仍然较低，表明脑组织与训练/主测试组织之间存在较强分布差异。
4. 这些结果是外部泛化测试结果，不能并入原七切片 LODO 主结果均值。

## 产物

- 评估脚本：`SCP-NET-work/experiments/external_generalization_strict_calibration_20261031/evaluate_external_strict_calibration.py`
- 逐数据集汇总：`SCP-NET-work/experiments/external_generalization_strict_calibration_20261031/external_summary.csv`
- 源验证选择网格：`SCP-NET-work/experiments/external_generalization_strict_calibration_20261031/source_selection_grid.csv`
- 预测结果：同一实验目录下的 `*_proteinconvex_strict_rna_predictions.npy`
