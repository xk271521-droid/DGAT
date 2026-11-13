# 三个跨组织目标指标审计

审计日期：2026-11-02

## 结论

截图中的 ProteinConvex 数值来自其 `results.json` 和逐蛋白指标文件。使用 ProteinConvex 每个方向 `results.json.common_panel.proteins` 指定的内部蛋白顺序，直接对保存的 prediction/truth 重新计算后，三个目标均能复现：

| 目标 | PCC-Pro | RMSE-Global |
|---|---:|---:|
| PBC_PR_6837 | 0.2438 | 0.2897 |
| PBC_PR_6835_5A | 0.2721 | 0.3621 |
| GBM | 0.2040 | 0.3176 |

因此，当前低分不是这三行指标在表格中被算错。ProteinConvex 在这三个组织外目标上的泛化确实明显弱于 NicheTrans 和 DPAS-Graph，尤其是两个 PBC 目标；GBM 上三个模型都较低，说明该组织的跨域难度更高。

## 已发现并修复的问题

ProteinConvex 预测列使用作者脚本根据第一个数据集变量顺序生成的共同蛋白顺序。此前整理 truth 文件时误用了 JSON 中的字母顺序。这样用错误 truth 文件独立复算，会得到不属于同一蛋白列的错误低分。

已按每个方向 `results.json.common_panel.proteins` 的真实顺序重新导出全部 18 个 ProteinConvex `target_source_selected_truth.npy`，并更新跨组织结果哈希和清单。预测文件、模型权重和正式 `results.json` 指标没有修改。

NicheTrans 和 DPAS-Graph 的三个目标已用各自保存的 prediction/truth 独立复算，和表格数值一致。

## 解释边界

ProteinConvex 的三组跨组织任务都遵守整组留出：PBC 任务 source 不含三个 Breast Cancer 数据集，GBM 任务 source 不含 GBM。当前结果反映的是模型跨组织泛化能力，而不是目标标签泄漏或表格列错位。
