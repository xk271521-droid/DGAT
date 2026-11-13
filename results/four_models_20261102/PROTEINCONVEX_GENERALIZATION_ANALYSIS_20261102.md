# ProteinConvex 泛化表现分析

日期：2026-11-02

## 已确认结果

ProteinConvex 的同组织与跨组织结果来自统一 9384/31 面板和 source-only checkpoint 选择。关键 PCC-Pro 如下：

| 实验 | 目标 | PCC-Pro |
|---|---|---:|
| 同组织 `A1_LN -> D1_LN` | D1_LN | 0.4772 |
| 同组织 `D1_LN -> A1_LN` | A1_LN | 0.5467 |
| 同组织 `PBC_PR_6835_5A -> PBC_PR_6837` | PBC_PR_6837 | 0.3323 |
| 同组织 `PBC_PR_6837 -> PBC_PR_6835_5A` | PBC_PR_6835_5A | 0.3492 |
| 跨组织留出 LN | A1_LN / D1_LN | 0.3317 / 0.3225 |
| 跨组织留出 Tonsil | 4 个目标 | 0.4704、0.4628、0.4811、0.6125 |
| 跨组织留出 Breast Cancer | Breast_Cancer | 0.5920 |
| 跨组织留出 Breast Cancer | PBC_PR_6837 / PBC_PR_6835_5A | 0.2438 / 0.2721 |
| 跨组织留出 GBM | GBM | 0.2040 |

## 解释

### D1_LN

`A1_LN -> D1_LN` 为 0.4772，反向 `D1_LN -> A1_LN` 为 0.5467。两条 source validation PCC 分别约为 0.649 和 0.658，说明训练侧验证表现接近，但目标切片表现有方向性差异。这更像切片间分布差异和方向不对称，不是输入维度或指标计算错误。

### PBC 两个目标

PBC 的问题在跨组织前已经出现：同组织双向只有 0.3323 和 0.3492。因此 PBC 的低分不能只归因于“跨组织留出”；它是 PBC 外部切片本身的域差异或测量差异。PBC 数据说明中明确记录它们是不同 slide，并且空间坐标尺度和当前 Breast section 不同。跨组织训练没有看到 PBC 的 ADT 标签，所以模型无法针对该域做校准。

### GBM

GBM 只有一个数据集，没有同组织配对，只能做整组留出。三个作者模型在 GBM 上的 PCC-Pro 也只有约 0.217、0.245 和 0.204，说明 GBM 是整个任务中较难的跨域目标，ProteinConvex 的 0.204 并不是单独异常。

### Breast_Cancer 与 PBC 的差异

同一个 Breast Cancer 留出组中，当前 `Breast_Cancer` 目标为 0.5920，但两个 PBC 目标为 0.2438 和 0.2721。这说明“留出 Breast Cancer”不是一个同质难度组：当前 Breast section 可以泛化，PBC 外部 sections 明显更难。表中按目标数据集分开评分是必要的，不能用一个 Breast Cancer 组均值掩盖 PBC 退化。

## 当前判断

目前没有证据表明 ProteinConvex 这几行是指标算错。已修复过一次 truth 列顺序问题后，ProteinConvex 的 prediction/truth 已按模型内部蛋白顺序保存，结果可以独立复现。现有证据支持“模型在 PBC 和 GBM 外部域上泛化不足”，而不是“表格计算错误”。

## 后续诊断顺序

1. 对 D1_LN、PBC 两个目标和 GBM 做逐蛋白 PCC、预测均值/方差与真实值均值/方差对比。
2. 对 source 与 target 的 RNA 共同基因分布做无标签 shift 检查，确认哪些基因驱动域差异。
3. 单独评估目标 RNA 的尺度、邻域图和 slice-context 描述是否落在 source 范围外。
4. 在不改变正式测试结果的前提下，做 source-only 的外部域校准或留一数据集诊断；未完成前不宣称改进模型在跨组织上优于作者模型。
