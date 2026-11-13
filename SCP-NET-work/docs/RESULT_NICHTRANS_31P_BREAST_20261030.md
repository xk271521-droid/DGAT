# NicheTrans-31P 独立架构适配：Breast_Cancer

本候选使用官方 NicheTrans 的 `model/nicheTrans.py` 核心结构：RNA MLP、中心/邻居语义 token、自注意力、GEGLU FFN 和逐蛋白预测头。将官方 Breast 示例的 313→2 适配为当前冻结面板的 9764→31，并重写为项目的 7 切片 LODO、source-only validation PCC-Pro 选点协议。

官方源码保持不变；当前适配使用 `noise_rate=0.2`、`dropout_rate=0.2`、40 epoch 和 NicheTrans 风格邻居腐蚀。

| 指标 | ProteinConvex 基线 | NicheTrans-31P |
|---|---:|---:|
| 源验证最佳 PCC-Pro | 0.751729 | 0.720202 |
| Breast PCC-Pro | 0.602035 | 0.548614 |
| Breast PCC-Spot | 0.633267 | 0.503515 |
| Breast RMSE | 0.299574 | 0.329438 |

该独立架构明显退化，不扩展七折。适配代码和运行目录已清理；官方 NicheTrans 源码仍保留在 `upstream/NicheTrans/`。
