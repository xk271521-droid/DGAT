# 清理记录

2026-11-02 完成三模型结果整理。

保留：

- ProteinConvex 最新十数据集正式运行目录；
- NicheTrans fixed-panel 十数据集正式运行目录；
- DPAS-Graph all-results 十数据集正式结果目录；
- 三模型比较工作簿、协议文档、官方源码和适配脚本。

删除：

- NicheTrans 初始运行、旧 rerun 目录；
- DPAS-Graph 初始运行、smoke test、v2 和服务器中间目录；
- ProteinConvex 旧七数据集运行、旧种子运行和单折旧目录；
- 旧四模型/协议比较表及其检查输出；
- 已完成上传后不再使用的 DPAS-Graph 传输压缩包和旧传输目录；
- 旧表格导出、测试、检查和临时传输脚本。

本次不删除 SPINE、DGAT 或其他历史研究线的结果和脚本。

整理后复核发现 ProteinConvex 的 8 个同组织目录缺少单独的 truth 文件。已从十数据集固定 31 蛋白面板重新导出 `target_source_selected_truth.npy`，并逐方向核对其形状与 prediction 文件一致；没有重新训练模型。
