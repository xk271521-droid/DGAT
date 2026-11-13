# 四篇论文的官方源码归档

下载日期：2026-11-01（America/Halifax）。本目录独立于SCP-NET和冻结数据。

| 项目 | 官方仓库 | 本地目录 | 固定提交 | 当前公开跟踪文件 |
|---|---|---|---|---:|
| DGAT | https://github.com/osmanbeyoglulab/DGAT | `DGAT/` | `b0128f342de70fa04546d27f9967eab1821f0933` | 58 |
| DPAS-Graph | https://github.com/xmy413/DPAS-Graph | `DPAS-Graph/` | `c3f80676051af5aa15969c9726c6d8bdf6db12f8` | 33 |
| NicheTrans | https://github.com/YSTLab/NicheTrans | `NicheTrans/` | `85fab4722f9b16cdef8b23f5a9ace0c15980ca16` | 53 |
| SPINE | https://github.com/VitaIntelli-CQU/SPINE | `SPINE/` | `961d765e817573222a18d0cc135d38aadce03c09` | 35 |

仓库来源分别由DGAT论文的Code availability及DPAS-Graph论文的Data availability明确给出：

- DGAT预印本：https://pmc.ncbi.nlm.nih.gov/articles/PMC12265741/ ，DOI `10.1101/2025.07.05.662121`。
- DGAT正式版：https://www.nature.com/articles/s41467-026-73114-z 。
- DPAS-Graph：https://doi.org/10.1093/bib/bbag480 。

四个项目均保留作者当前公开源码，没有修改参数、目录结构、预处理代码或训练代码。
每个跟踪文件都与Git对象逐字节校验；工作目录干净，`git fsck --full`通过。
没有未取回的Git子模块或LFS指针。详细版本、每文件SHA256及归档校验值见 `SOURCE_MANIFEST.json` 及各项目清单。

“完整”指上述固定提交的全部公开跟踪文件，不表示已获得作者未公开文件、外部数据、
预训练权重或论文每个实验的全部运行产物，也不表示已验证能够直接复现论文分数。

## 版本和离线保存

- DGAT：当前版本完整源码，Git历史为浅克隆，只保存当前提交；不是缺少当前源码文件。
- DPAS-Graph：完整普通Git克隆，非浅克隆。
- `archives/DGAT_b0128f342de70fa04546d27f9967eab1821f0933.zip`：GitHub codeload官方固定提交压缩包。
- `archives/DPAS-Graph_c3f80676051af5aa15969c9726c6d8bdf6db12f8.zip`：从已校验官方Git树生成的源码包。
- `archives/DGAT_verified_archive_extract_b0128f3/`：曾用于逐字节对照的官方压缩包解压副本。
- `download_attempts/`：中断的完整历史下载元数据，不是可用源码目录。

后续服务器不能访问外网时，可以传源码压缩包和本目录清单，不需要重新从非官方镜像找代码。

## 常用入口

- DGAT训练示例：`DGAT/Demo1_Train.ipynb`；多样本训练：`DGAT/Pretrain_DGAT.ipynb`。
- DGAT网络：`DGAT/Model/`；论文图复现：`DGAT/Reproduction/`。
- DPAS-Graph留一训练：`DPAS-Graph/scripts/train/run_lodo.py`。
- DPAS-Graph网络：`DPAS-Graph/dpas-graph/models/dpas_graph.py`。
- DPAS-Graph训练实现：`DPAS-Graph/dpas-graph/training/trainer.py`。
- NicheTrans官方代码：`NicheTrans/`；主要模型入口：`NicheTrans/model/`，乳腺癌训练工具：
  `NicheTrans/utils/utils_training_breast_cancer.py`，参数文件：`NicheTrans/args/args_breast_cancer.py`。
- SPINE训练入口：`SPINE/spine/app/flow/train_rna_to_protein.py`；模型实现：`SPINE/spine/model/`；
  数据与预处理：`SPINE/spine/data/`、`SPINE/spine/app/preprocessing/`。

本次只下载和静态检查，未安装作者依赖、未执行作者脚本。
DPAS-Graph公开目录的导入路径问题及预处理差异见 `UPSTREAM_USAGE_NOTES.md`，没有偷偷修补作者文件。

## 外部文件

DGAT README把训练数据、预训练权重和论文复现结果放在作者Google Drive：
https://drive.google.com/drive/folders/1M9gIOFjK4wBk7DEj9MrwAlV8OHEgS2BA 。
这些外部数据和权重不包含在Git仓库，本次没有额外下载。
DPAS-Graph README说明原始数据从10x和GEO获取，没有随源码分发。
当前项目已冻结的数据继续原样保留，本次源码下载不触碰它。

NicheTrans官方仓库固定提交为 `85fab4722f9b16cdef8b23f5a9ace0c15980ca16`，
已原样归档53个公开文件（约9.65 MB），逐文件校验见 `NicheTrans_SOURCE_MANIFEST.json`。
仓库README指向的处理后数据、外部数据和预训练结果没有下载；本次只保存代码和notebook。

SPINE官方仓库固定提交为 `961d765e817573222a18d0cc135d38aadce03c09`，
已原样归档35个公开文件（约0.87 MB），逐文件校验见 `SPINE_SOURCE_MANIFEST.json`。
仓库README说明的demo数据未随源码包下载；本次只保存官方源码和依赖说明。

## 下载过程

直连GitHub发生连接重置或超时，随后通过电脑现有代理访问官方GitHub；没有使用未知第三方代码包，
也没有修改系统或Git全局代理。DGAT完整历史拉取中断后改取当前版本，浅克隆最终正常完成，
且与固定提交官方压缩包的58个文件全部逐字节相同。
DPAS-Graph普通克隆正常完成。两仓库本地Git配置将autocrlf设为false以保留原始文件字节。
NicheTrans因Git clone和codeload连接不稳定，改用官方固定提交的raw文件逐个下载并按Git blob SHA校验。
SPINE直连Git clone发生连接重置，改用官方固定提交的GitHub codeload压缩包，并按固定提交重新解压。
未安装任何作者依赖、未执行作者脚本、未连接实验服务器。
