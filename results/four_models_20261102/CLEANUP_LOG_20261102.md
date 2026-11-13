# 旧乳腺癌泄漏结果清理记录

日期：2026-11-02

按用户要求删除 NicheTrans 和 DPAS-Graph 错误 source 划分产生的旧 Breast Cancer 留出结果：

- 本地 `SCP-NET-work/runs/nichetrans_10_dataset_20261101_fixed_panel/cross_holdout_Breast_Cancer/`
- 本地 `SCP-NET-work/runs/dpas_graph_10_dataset_20261101_all_results/cross_holdout_Breast_Cancer/`
- 服务器 `/public/liqiang/DGAT/SCP-NET-work/runs/dpas_graph_10_dataset_20261101_server/cross_holdout_Breast_Cancer/`

四模型工作簿里对应 NicheTrans/DPAS-Graph 的 6 行已清空旧数值，并标记为 `已删除 / 待重跑`。桌面副本 `C:\Users\xk2714\Desktop\four_model_comparison_20261102.xlsx` 已替换为同一份更新后工作簿；它与项目输出副本 SHA256 一致。

保留的修正版：

- NicheTrans：`SCP-NET-work/runs/nichetrans_10_dataset_20261102_corrected_breast_holdout/cross_holdout_Breast_Cancer/`
- DPAS-Graph：服务器 `/public/liqiang/DGAT/SCP-NET-work/runs/dpas_graph_10_dataset_20261102_corrected_breast_holdout/cross_holdout_Breast_Cancer/`

DPAS-Graph 修正版后台训练进程仍在运行；本次退出监控不停止该进程。
