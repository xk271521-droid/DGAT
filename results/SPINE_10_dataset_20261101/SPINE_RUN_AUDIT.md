# SPINE 运行审计记录（2026-11-02）

## 可以确认的训练配置

从服务器保存的每个 `config.json` 核对到的配置一致：

- epochs 上限 500，early_stop_patience 50，early_stop_metric=pcc
- batch_size=1，gradient_accumulation_steps=2
- lr=5e-5，weight_decay=1e-4，clip_norm=1.0
- RNA feature_dim=9384，protein 数=31
- hidden_dim=256，pairwise_hidden_dim=256，n_layers=6，n_heads=8
- dropout=0.3，attn_dropout=0.3
- n_neighbors=8，MLP intermediate=1024，MLP layers=2
- sample_times=10，n_sample_steps=10，prior_sampler=gaussian
- gene_recon_weight=0.2，loss_type=mse，seed=42
- normalize_method=None，rna_feature_normalize=None

这些值与本地 `upstream/SPINE/spine/app/flow/train_rna_to_protein.py` 当前入口中的默认值/固定值一致；没有证据表明训练时把 epoch、学习率、模型层数或早停耐心改小。

## 本次适配实际改变的内容

1. 没有直接修改 `upstream/SPINE` 的模型和训练源码；`train_rna_to_protein_adapted.py` 只增加 NumPy 标量 JSON 序列化，并把输出目录从带时间戳改为固定实验目录，便于后台队列管理。
2. `prepare_spine_direction.py` 是我们新增的数据适配器：从统一的 10 数据集 h5ad 中取跨所有数据集共有的 9384 个 RNA 基因和 31 个蛋白，生成 SPINE 所需的 RNA_EMBED、protein.h5ad 和 split 文件。
3. 每个 source 组织内部按 seed=42 随机抽 10% 做 validation；target 组织不参与 checkpoint 选择。这个 split 是本项目实验协议，不是作者 demo 的原始 split 文件。
4. 因此，当前结果可以描述为“SPINE 模型在本项目统一输入/同组织与跨组织协议上的复现运行”，不能直接声称是作者原始数据预处理和原始 split 的逐字复现。

## 为什么比“两天”快

服务器实际用了两张 GPU 并行跑任务；单个方向达到早停通常约 135–193 epoch，两个 Human Tonsil 方向跑到 500 epoch。任务按队列并行后，总墙钟时间会明显短于单卡逐个跑满 500 epoch 的估计。师兄说的“两天”应是资源、并行度、数据规模或保守估计不同；它不是 SPINE 必须固定运行两天。

## 结果使用边界

最终预测文件已经保存并拉回本地。服务器 v5 的 Breast Cancer 跨组织方向曾因错误 source 划分缺少两个真正的 held-out 输入；该旧结果已废弃。服务器 v6 已按乳腺癌组织整体留出重新训练，并生成 `Breast_Cancer`、`PBC_PR_6837`、`PBC_PR_6835_5A` 三个目标的完整预测。v6 审计材料保存在 `results/SPINE_10_dataset_20261101/audit/cross_holdout_Breast_Cancer_v6/`。

由于本地 SPINE 副本没有可用的 git remote/commit 元数据，当前无法仅凭本地文件证明它对应作者仓库的某个具体 commit。若论文需要严格作者版本复现，应补记作者仓库 commit 或重新从官方发布版本锁定后跑一套审计实验。
