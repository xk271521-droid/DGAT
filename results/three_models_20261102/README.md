# 三模型正式结果

整理日期：2026-11-02

状态：该三模型工作簿已于 2026-11-03 被 `results/five_models_20261103/comparison/five_model_comparison_20261103.xlsx` 取代。原始三模型运行结果目录仍是四模型/五模型表共用的数据来源。

本目录对应当前十数据集正式比较中的三个已完成模型：

- `ProteinConvex + Strict RNA`
- `NicheTrans`
- `DPAS-Graph`

统一协议：10 个数据集、9384 个共同 RNA 基因、31 个共同 ADT 蛋白；同组织 8 个双向方向，跨组织 4 个整组留出方向。每个模型均应有 12 个方向、18 个目标结果。

正式运行结果仍保存在各自的运行目录中，避免复制大体积 checkpoint 和中间缓存：

- ProteinConvex：`SCP-NET-work/runs/protein_convex_10_dataset_20261101`
- NicheTrans：`SCP-NET-work/runs/nichetrans_10_dataset_20261101_fixed_panel`
- DPAS-Graph：`SCP-NET-work/runs/dpas_graph_10_dataset_20261101_all_results`

本目录的 `RESULTS_MANIFEST.json` 记录结果覆盖和来源位置；`comparison/three_model_comparison_20261102.xlsx` 是三个模型的统一比较表。

跨组织低分及 truth 列顺序问题的审计记录见 `METRIC_AUDIT_20261102.md`。

每个正式运行目录保留最终预测、真实值、指标、`results.json`、配置和必要审计文件。旧版运行、smoke test、重复导出表和临时传输包已清理，SPINE、DGAT 及其他历史研究线不在本次清理范围内。
