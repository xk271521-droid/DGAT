# CURRENT_BENCHMARK (10 数据集 18 目标最终基准总览)

更新日期：2026-11-05
当前官方主导模型：`ProteinConvex + 初始粗糙度尺度门控 + Strict RNA + Source Affine`

---

## 一、 六模型全量 18 目标等权宏平均总览

| 模型类别 | 模型全称 | 18 目标宏平均 PCC-Pro | 18 目标全局 RMSE-Global | 跨组织 10 目标 PCC | 跨组织 10 目标 RMSE |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **Ours (最佳版本)** | **ProteinConvex + 初始粗糙度尺度门控 + Strict RNA + Source Affine** | **0.4695** | **0.2886** | **0.4178** | **0.3116** |
| **Ours (消融基线)** | ProteinConvex + Strict RNA Calibration (原版) | 0.4634 | 0.3073 | 0.4068 | 0.3344 |
| **Competitor 1** | NicheTrans | 0.4408 | 0.3082 | 0.3706 | 0.3212 |
| **Competitor 2** | DPAS-Graph | 0.4300 | 0.3204 | 0.3768 | 0.3393 |
| **Competitor 3** | SPINE | 0.4089 | 0.3175 | 0.3632 | 0.3342 |
| **Competitor 4** | DGAT | 0.4150 | 0.3262 | 0.3643 | 0.3461 |

---

## 二、 关键产物对应关系
- **正式表格**: `C:\Users\xk2714\Desktop\six_model_comparison_20261104_retained_clean.xlsx`
- **模型权重**: `SCP-NET-work/runs/relative_roughness_gated_runs/`
- **校准预测**: `SCP-NET-work/runs/relative_roughness_gated_runs_affine_calibrated/`

> 复现说明：本表指标对应当前归档权重的实际执行脚本。same 任务使用 `alpha=0`；cross 任务在 `initial_roughness < 0.008` 时使用 `alpha=15` 并惩罚当前绝对粗糙度，否则 `alpha=0`。正式模型名称按该实际逻辑统一为“初始粗糙度尺度门控”。总览 MAE-Global 已按 18 个方向明细更正为 `0.1852496849`。详见 `REPRODUCIBILITY_AUDIT_20261104.md`。
