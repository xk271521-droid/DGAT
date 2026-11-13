# 服务器上传大小清单

日期：2026-11-01

## 总体大小

| 目录 | 大小 | 用途 |
|---|---:|---|
| `E:\DGAT\SCP-NET-work` | 792.5 MiB | 项目代码、实验脚本、文档和历史/当前运行结果 |
| `E:\DGAT\ST2P` | 2.42 GiB | 数据、预处理矩阵、源码和实验归档 |
| `E:\DGAT\upstream` | 78.3 MiB | 四篇论文官方源码和源码归档 |
| 整个 `E:\DGAT` | 约 3.3 GiB | 当前工作区总量 |

## 服务器运行 ProteinConvex 所需

当前 10 数据集 ProteinConvex 训练只需要：

- 10 个数据集的独立预处理矩阵：约 1.52 GiB；
- `ST2P/src/`：约 0.3 MiB；
- `SCP-NET-work/experiments/protein_convex/` 和 `protein_convex_10_dataset_20261101/`：约数 MiB；
- 当前训练运行结果和日志：随运行增加，目前约数百 MiB以内。

因此，服务器如果只继续跑当前 ProteinConvex，建议先上传约 **1.6 GiB** 的“预处理数据 + 代码 + 运行配置”包，不需要上传全部历史结果。

## 后续五模型统一实验所需

如果服务器还要运行 DGAT、DPAS-Graph、NicheTrans、SPINE，额外上传 `upstream/`，约 **78.3 MiB**。但四个作者环境不能直接共用当前 Windows `st2p` 环境，应在服务器按各自 requirements/environment 文件建立独立环境。

## 暂时不必上传

- `ST2P/data/raw/` 原始矩阵：当前训练使用已完成的独立预处理矩阵；
- `SCP-NET-work/runs/` 中与新 10 数据集无关的旧运行结果；
- `__MACOSX/`、`tmp_spreadsheet_edit/` 和临时文件；
- 旧的 9764 基因外部映射预测结果：已废止。

## 注意

当前 10 数据集目录使用 Windows junction 指向已有数据目录。上传时不要只压缩 junction 本身，必须把实际目标目录中的文件复制到服务器，或者按本清单分别上传：

- `ST2P/data/processed/DGAT_final7_breast_ptprc31_muon017_axis0_20261022/`
- `ST2P/data/processed/PBC_breast_external_20261030/`
- `ST2P/data/processed/GBM_external_20261030/`
- `ST2P/data/processed/DGAT_final10_common9384_muon017_axis0_20261101/prepared_manifest.json`
