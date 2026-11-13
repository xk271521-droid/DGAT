# 官方源码使用记录（未改作者文件）

日期：2026-10-30。仅静态阅读，没有执行作者程序。

## DPAS-Graph目录与导入名不一致

固定提交：`c3f80676051af5aa15969c9726c6d8bdf6db12f8`。

- 仓库包含 `dpas-graph/`，没有 `src/dpas/`。
- `scripts/train/run_lodo.py:5` 和 `scripts/eval/benchmark_inference.py:5` 将 `repo/src` 加入搜索路径。
- 训练入口使用 `from dpas.training.trainer import train_and_evaluate_fold`。
- 模型和数据模块内部也使用 `dpas.*` 绝对导入。
- 当前跟踪文件中没有setup.py、pyproject.toml或另一个dpas包来解释这个映射。

因此直接照README运行存在找不到dpas包的明确静态风险。这是公开仓库原有布局，
本次完整下载不等于已解决运行问题。后续需要在独立工作副本做路径适配并保留差异记录；
不在官方归档里直接改名或改源码。

## 预处理不能直接当成与本项目完全一样

`DPAS-Graph/scripts/train/run_lodo.py` 的 `_adt_clr_inplace` 在raw分支执行
`log1p(X)` 后减去每个spot跨蛋白的均值（axis=1）。若输入标记为proc则跳过此变换。
输入文件名含proc或clr时，入口会把它标记为proc。

当前SCP-NET冻结注册表使用Muon0.1.7 CLR axis=0、保留其非负值尺度。
这两个实现不能不经核对就称为同一个变换。DPAS-Graph仓库还有独立的data_prep脚本，
需要连同输入文件状态一起判断；本记录没有认定论文实际每个实验都走了哪个分支。
后续比较必须明确使用哪套冻结输入，不能为追分静默重算目标数据。

本次没有改动任何一套预处理，也没有把DPAS论文表格数值和本项目结果直接当作同尺度比较。

## DGAT外部运行材料

固定提交：`b0128f342de70fa04546d27f9967eab1821f0933`。
公开源码的训练示例、网络、工具和复现notebooks都已保留。
README指向Google Drive的数据、权重和论文复现结果属于额外材料，本次源码下载不包含它们。
当前仓库README写Python3.11；论文不同版本的软件环境说明可能不同，后续要按固定版本单独核对。

## 当前工作边界

本次下载用于后续代码审查和可追溯复现，不自动改变现有ProteinConvex的训练参数、
7切片面板、邻居构造、数据集或source-only checkpoint选择规则。

## DPAS-Graph公开训练入口的选点规则

2026-10-30继续静态检查发现：固定提交的
`dpas-graph/training/trainer.py` 约489行在每个epoch调用`evaluate_testset`，
约505行根据`test_metrics["rmse_global"]`更新`best_state`，最终加载该状态。
此代码路径使用目标集指标选择checkpoint，与本项目source-only验证规则不同。
这是一条特定公开代码路径的可观察行为，不据此猜测作者所有运行的未公开细节。
本项目借鉴RNA重建时仍保留source-validation PCC-Pro选点，不能通过改成目标选点来宣称改进。

## NicheTrans官方源码归档

固定仓库：`https://github.com/YSTLab/NicheTrans`。
固定提交：`85fab4722f9b16cdef8b23f5a9ace0c15980ca16`。
本地目录：`upstream/NicheTrans/`；文件清单和逐文件SHA256/Git blob校验见
`upstream/NicheTrans_SOURCE_MANIFEST.json`。

该提交包含53个公开文件，约9.65 MB，主要为模型、数据管理器、训练工具、参数文件和notebook。
README指向的外部数据和预训练材料不在代码仓库，本次没有下载。仅完成源码归档和静态阅读，
没有安装依赖、运行训练或把NicheTrans参数替换到当前SCP-NET/ProteinConvex协议中。

## SPINE官方源码归档

官方仓库：`https://github.com/VitaIntelli-CQU/SPINE`。
固定提交：`961d765e817573222a18d0cc135d38aadce03c09`。
本地目录：`upstream/SPINE/`；文件清单和归档哈希见 `upstream/SPINE_SOURCE_MANIFEST.json`。

该提交包含35个公开文件，包含训练入口、flow、denoiser、transformer、数据加载和预处理代码。
README推荐Python 3.12.9和PyTorch 2.7.0+cu118；本次没有安装作者依赖、执行作者训练或下载demo数据。
官方源码保持原样，后续如需适配10数据集，只在独立适配层中处理并保留差异记录。
