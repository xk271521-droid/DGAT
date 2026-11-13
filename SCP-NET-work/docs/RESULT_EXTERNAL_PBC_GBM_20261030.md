# 外部切片验证：两个新乳腺癌病人 + 胶质瘤（2026-10-30）

> 历史记录：本结果使用旧的外部映射输入。新的 10 数据集实验必须重新使用 9384 个共同 RNA 基因面板。

目的：检验现有 Breast 折 0.602 是否能跨病人、跨组织泛化。外部切片只用 RNA + 空间坐标，目标 ADT 只在预测完成后读取一次。

## 数据与协议

- 外部切片：
  - `PBC_PR_6837`：4,850 spots（DGAT 官方训练数据，与现有 Breast 不同切片）。
  - `PBC_PR_6835_5A`：3,243 spots。
  - `GBM`：5,637 spots（10x Human Glioblastoma，同 31 抗体面板）。
- 预处理：与冻结 7 切片完全相同的 DGAT QC / normalize / log1p / scale / clip 10，ADT 用 Muon-0.1.7 CLR axis=0。
- 基因面板：9,764。QC 后外部切片分别保留 9,477 / 9,595 / 9,594 个，缺失基因在推理时置 0（标准化空间均值）。
- 模型：
  - **zero-shot**：Breast 折 checkpoint（训练时完全没见过乳腺、PBC、GBM）。
  - **one-shot**：A1_LN 折 checkpoint（训练时见过 1 个 Breast 病人，没见过 PBC/GBM）。
- 管线自检：同一推理脚本在冻结 Breast 上复现 `0.6020347 / 0.6332673 / 0.2995744`，与正式结果一致。

## Zero-shot：4 个模型 × 3 个外部切片

| 切片 | baseline | CognateAnchor | ProteinConvex | MultiView | RNA-only 基线 |
|---|---:|---:|---:|---:|---:|
| PBC_PR_6837 | 0.2135 | **0.2337** | 0.2257 | 0.2256 | 0.1696 |
| PBC_PR_6835_5A | 0.2479 | **0.2662** | 0.2468 | 0.2381 | 0.1439 |
| GBM | 0.0973 | 0.1107 | **0.1262** | 0.1091 | 0.1369 |

Breast 同 checkpoint 参考：ProteinConvex `0.6020 / 0.6333 / 0.2996`。

## One-shot（训练时见过 1 个乳腺病人）

| 切片 | zero-shot | one-shot | 变化 | RNA 基线 |
|---|---:|---:|---:|---:|
| PBC_PR_6837 | 0.2257 | 0.2633 | +0.0376 | 0.1696 |
| PBC_PR_6835_5A | 0.2468 | 0.2902 | +0.0434 | 0.1439 |
| GBM | 0.1262 | 0.1609 | +0.0347 | 0.1369 |

## 结论

1. **0.602 不能跨病人复制。** 换成另外两个乳腺癌病人，绝对 PCC-Pro 只有 0.21–0.29；换成脑组织只有 0.10–0.16。单看一个 Breast 切片确实证明不了泛化性。
2. **模型不是完全没用。** 在 PBC 上仍高于 RNA-only 基线 +0.04–0.12；GBM zero-shot 低于 RNA 基线（−0.011），one-shot 才略高于基线（+0.024）。
3. **加一个乳腺训练病人只带来约 +0.04。** 这能改善外部病人，但远不足以把 0.23–0.29 拉到 0.60，说明外部切片与当前 Breast 的差距不只是“没见过这个病人”。
4. **外部切片本身更难/分布不同。** RNA-蛋白逐蛋白相关均值从 Breast 的 0.398 降到 PBC 的 0.17/0.14、GBM 的 0.137；绝对分数不能直接和 Breast 横比。
5. **模型排名在外部切片上不稳定。** PBC 上 CognateAnchor 最好，GBM 上 ProteinConvex 最好，PBC_6835 上 baseline 与 CognateAnchor 接近；“ProteinConvex 最好”只在当前 Breast 切片成立。

## 限制

- 外部切片 QC 后缺失 170–287 个 panel 基因，推理时置 0；不足以解释 0.60→0.23 的差距，但会轻微压低分数。
- PBC 与 GBM 的队列/建库/组织质量与当前 Breast 不同；低 RNA-蛋白相关说明数据本身更难。
- 每个模型只有一个 seed；外部切片的 seed 方差未测。
- PBC/GBM 是外部测试，不参与任何训练或选点。

## 决策

- 不再把 `0.602` 当作泛化性能；论文中应标注为“同队列/单切片”结果。
- 泛化证据改为：7 折跨组织 + 两个外部乳腺癌病人 + 脑组织外部测试，并同时报告 RNA-only 基线。
- 若目标是提高外部病人性能，必须把 PBC 加入训练并让所有模型在同一扩展源域重跑；按 one-shot 的结果，预期只有约 +0.04，不能承诺回到 0.60。

## 产物

- 预处理：`ST2P/data/processed/PBC_breast_external_20261030`、`ST2P/data/processed/GBM_external_20261030`。
- 推理脚本：`SCP-NET-work/experiments/pbc_external_validation_20261030/external_inference.py`。
- 结果：`SCP-NET-work/experiments/pbc_external_validation_20261030/results/`。
