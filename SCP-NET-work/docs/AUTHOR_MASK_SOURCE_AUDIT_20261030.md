# DGAT与DPAS-Graph官方代码的mask/dropout来源核对

核对日期：2026-10-30。核对对象是本地归档的两个官方代码目录：
`upstream/DGAT/` 和 `upstream/DPAS-Graph/`。

## 结论

在这两份本地官方代码中，没有找到“50% RNA feature mask/dropout”设置。
因此当前ProteinConvex运行参数`gene_dropout_rate=0.5`不能称为DGAT或DPAS-Graph作者的原始设置。

## DGAT

官方 `upstream/DGAT/Model/Train_and_Predict.py` 使用全局`dropout_rate=0.3`，
传给GAT encoder；`upstream/DGAT/Model/dgat.py`中对应的是隐藏层/图注意力特征上的普通`nn.Dropout`。
官方训练脚本的`epochs=100`，没有`gene_dropout_rate`、随机RNA列mask或50% spot mask。
官方脚本还在训练过程中每轮评估test并记录test指标，不能把该选点行为复制到当前严格的source-validation协议。

## DPAS-Graph

官方 `upstream/DPAS-Graph/dpas-graph/training/trainer.py` 同样使用`dropout_rate=0.3`。
`Decoder_mRNA_MaskedSmall`中的“Masked”是类名和一个可调用的`masked_mse`工具，
但当前官方trainer实际使用的是普通`rmse_loss(rna_recon, m_x)`，没有创建随机mask。

官方 `dpas-graph/data/graph_dataset.py` 有可选的边权dropout，但环境变量
`DPAS_EDGE_DROP_BASE`默认值为`0.0`，不是50%；它作用于图边权，不是RNA基因列或spot。

## 当前项目实际实现

ProteinConvex运行配置使用`gene_dropout_rate=0.5`和`neighbor_mask_rate=0.25`。
前者在 `ST2P/src/v9_robust_source_subset_test_oracle.py` 的`feature_dropout`中通过
`F.dropout(x, p=rate, training=True)`对中心、local和regional RNA输入逐元素随机置零；
后者随机屏蔽邻居token，并保证每个spot至少保留一个邻居。
这两项是当前项目的训练增强，不是已核实的DGAT或DPAS-Graph官方50%设置。
当前代码注释明确写的是“`.5 matches NicheTrans's noise dropout`”。

## 对当前分数问题的影响

这不能直接证明50% dropout导致了ProteinConvex分数低，也不能证明去掉它一定会提高分数。
但它说明之前把该参数解释成“作者为Breast泛化设计的50% mask”是不准确的。
在继续改架构前，应先做一个来源清楚的受控对照：保持ProteinConvex结构、数据、source-only选点和20 epoch不变，
单独比较当前额外RNA feature dropout/邻居token mask与官方代码实际使用的普通模型dropout设置。
不能把这个对照称为DGAT或DPAS-Graph论文的完整复现。

## 本地证据

- DGAT训练入口：`upstream/DGAT/Model/Train_and_Predict.py:22,87,101`。
- DGAT GAT dropout：`upstream/DGAT/Model/dgat.py:8,47,66,73`。
- DPAS训练入口：`upstream/DPAS-Graph/dpas-graph/training/trainer.py:25,401`。
- DPAS可选边权dropout默认值：`upstream/DPAS-Graph/dpas-graph/data/graph_dataset.py:19,154`。
- DPAS RNA decoder及实际普通RMSE：`upstream/DPAS-Graph/dpas-graph/models/dpas_graph.py:211,230`、
  `upstream/DPAS-Graph/dpas-graph/training/trainer.py:485`。
- 当前项目feature dropout：`ST2P/src/v9_robust_source_subset_test_oracle.py:120,205`。
