# DGAT / SCP-NET 最终项目归档索引 (PROJECT_INDEX)

更新日期：2026-11-05
实验基准：10 个完整空间转录组数据集，共同 RNA 9384 基因，共同 ADT 31 蛋白。

---

## 一、 核心保留模型轨道 (Retained Models)
本项目严格锁定并保留以下核心模型与对照基准数据，所有产物与预测矩阵均完好留存：

1. **全新最佳自适应模型 (Ours - Final Best)**:
   - **名称**: `ProteinConvex + 初始粗糙度尺度门控 + Strict RNA + Source Affine`
   - **核心训练代码**: `SCP-NET-work/experiments/relative_roughness_gated_selection/scn_pro_relative_roughness_gated.py`
   - **校准脚本**: `SCP-NET-work/experiments/relative_roughness_gated_selection/run_source_only_rna_affine_calibration.py`
   - **训练 Checkpoint**: `SCP-NET-work/runs/relative_roughness_gated_runs/` (包含全部 18 个测试目标的 best_source_validation.pth)
   - **最终校准预测产物**: `SCP-NET-work/runs/relative_roughness_gated_runs_affine_calibrated/`
   - **核心表现**:
     - 18 目标宏平均 PCC-Pro: **0.4695** (全场第一)
     - 18 目标全局 RMSE-Global: **0.2886** (全场第一)
     - 18 目标宏平均 MAE-Global: **0.1852**
     - 跨组织 10 目标均 PCC-Pro: **0.4178**, RMSE-Global: **0.3116** (8/10 目标同时战胜全部对手)
     - 同组织 8 目标均 PCC-Pro: **0.5341**, RMSE-Global: **0.2598** (8/8 目标全部斩获第一)

2. **基线消融模型 (Ablation Baseline)**:
   - **名称**: `ProteinConvex + Strict RNA Calibration`
   - **路径**: `SCP-NET-work/runs/protein_convex_rna_calibration_10_dataset_20261103/`
   - **作用**: 用于论文消融实验，证明未加自适应门控与仿射校准前的原始效果。

3. **四个标准公开对比基线 (Competitors)**:
   - **NicheTrans**: `SCP-NET-work/runs/nichetrans_10_dataset_20261102_corrected_breast_holdout/`
   - **DPAS-Graph**: `SCP-NET-work/runs/dpas_graph_10_dataset_20261101_all_results/`
   - **SPINE**: `SCP-NET-work/runs/protein_convex_10_dataset_20261101/` (同面板 SPINE 对照指标)
   - **DGAT**: `SCP-NET-work/runs/dgat_10_dataset_20261101_all_results/`

---

## 二、 关键对比大表
- **桌面最新干净六模型大表**: 
  `C:\Users\xk2714\Desktop\six_model_comparison_20261104_retained_clean.xlsx`
  (已删除未校准的原始冗余行，只保留你圈定的 5 个标准模型 + 我们的终极改进模型，表格干净整洁，格式完全规范)

---

## 三、 清理与归档记录
已彻底清除以下历史冗余/失败的试跑目录与代码：
- 清理一刀切粗糙度失败产物 (`auto_roughness_strict_calibration_10_dataset_20261103`)
- 清理复合域距离选择器试跑产物 (`source_only_affine_selector_v3 / v4`)
- 清理历史平滑临时脚本与无用 dump 文件
工作区已保持绝对纯净与高可用状态。

## 四、复现口径审计

版本名称按“ProteinConvex + 初始粗糙度尺度门控 + Strict RNA + Source Affine”锁定。现有 18 个方向权重和指标对应的实际训练逻辑为：same 方向 `alpha=0`；cross 方向在 `initial_roughness < 0.008` 时使用 `alpha=15`，否则为 0，并惩罚当前绝对粗糙度。`relative_roughness_gated_selection` 是历史目录名，不代表本批产物使用了相对增长率公式。桌面正式表的总览 MAE-Global 已按 18 个方向明细更正为 **0.1852496849**。详见 `SCP-NET-work/docs/REPRODUCIBILITY_AUDIT_20261104.md`。

## 五、私有源码仓库（2026-11-12）

- 用户自己的私有仓库：`https://github.com/xk271521-droid/DGAT`。
- 本地发布副本：`E:/DGAT/git-publication/`；不改变旧 `ST2P/.git` 的远程或历史。
- 包含当前改进模型、共享依赖、其他保留源码、运行源码快照、实验配置、紧凑指标和文档，
  以及四个作者固定提交的全部公开跟踪文件。
- 数据矩阵、权重、预测/注意力数组、环境和临时文件留在本地。
- 根 `README.md` 记录代码入口和运行限制；
  `SCP-NET-work/docs/GIT_SOURCE_ARCHIVE_20261112.md` 记录归档口径；
  发布目录 `GIT_SOURCE_MANIFEST.json` 保存逐文件 SHA-256。
- 推送状态和提交哈希见根目录 `GIT_UPLOAD_RECEIPT_20261112.json`。
