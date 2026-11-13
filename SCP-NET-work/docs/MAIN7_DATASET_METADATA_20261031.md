# 7 个主数据集组织与平台元信息表

日期：2026-10-31

## 结论先行

本地原始注册表和 ST2P 官方仓库的描述都是：

```text
2 LN + 4 Tonsil + 1 Breast Cancer
```

因此，7 个主数据集不能简单按“A1/D1 四个 vs 剩下三个”直接称为同组织分组。正式分组实验前，还要确认每个切片的具体平台、病人/样本和批次信息。

官方证据：

- ST2P README 的冻结 7-section benchmark 说明；
- `src/build_final7_breast_raw_registry.py` 中的 registry 字段明确写为 `2 LN, 4 Tonsil, 1 Breast Cancer`；
- `docs/REPRODUCIBILITY.md` 列出了 7 个切片名称。

## 数据集清单

| 数据集 | 组织类别（本地注册表） | 样本/病人信息 | 平台/检测信息 | 本地已确认事实 | 同组织分组状态 |
|---|---|---|---|---|---|
| `A1_LN` | Lymph node（LN） | 未在当前本地 manifest 中记录 | 精确平台待核实；raw H5AD 含 `gene_ids/feature_types/genome` | 3,484 raw spots；3,293 QC 后；31 ADT | 与 `D1_LN` 属同一 LN 类别；平台/病人待核实 |
| `D1_LN` | Lymph node（LN） | 未在当前本地 manifest 中记录 | 精确平台待核实；raw H5AD 含 `gene_ids/feature_types/genome` | 3,359 raw spots；2,487 QC 后；31 ADT | 与 `A1_LN` 属同一 LN 类别；平台/病人待核实 |
| `A1_TNSL` | Tonsil（TNSL） | 未在当前本地 manifest 中记录 | 精确平台待核实；raw H5AD 使用 `feature_id/feature_type` | 4,337 raw spots；4,323 QC 后；31 ADT | 与其他 3 个 Tonsil 类数据集同组织类别；平台/病人待核实 |
| `D1_TNSL` | Tonsil（TNSL） | 未在当前本地 manifest 中记录 | 精确平台待核实；raw H5AD 使用 `feature_id/feature_type` | 4,519 raw spots；4,482 QC 后；31 ADT | 与其他 3 个 Tonsil 类数据集同组织类别；平台/病人待核实 |
| `HumanTonsil_10x` | Tonsil | 未在当前本地 manifest 中记录 | 数据集名称明确含 `10x`；完整平台细节待核实 | 4,194 raw spots；4,191 QC 后；31 ADT | Tonsil 类；需要确认是否与 TNSL/AddOns 可直接合组 |
| `HumanTonsil_AddOns_S1` | Tonsil | 未在当前本地 manifest 中记录 | 数据集名称明确含 `AddOns_S1`；完整平台细节待核实 | 4,313 raw spots；4,311 QC 后；原始 41 ADT，统一后 31 ADT | Tonsil 类；需要确认是否与 TNSL/10x 可直接合组 |
| `Breast_Cancer` | Breast cancer | 未在当前本地 manifest 中记录 | raw registry 明确来自 10x CytAssist FFPE Protein Expression 文件 | 4,169 raw spots；4,160 QC 后；31 ADT | 独立 Breast 类，不能与 LN/Tonsil 直接称同组织 |

## 当前可以确认的组织分组

### LN 类别

- `A1_LN`
- `D1_LN`

当前只能确认它们属于同一 LN 组织类别。是否同病人、同平台、同批次，现有 manifest 没有记录。

### Tonsil 类别

- `A1_TNSL`
- `D1_TNSL`
- `HumanTonsil_10x`
- `HumanTonsil_AddOns_S1`

当前注册表将它们归入 Tonsil 类别，但 `10x`、`AddOns_S1` 和 TNSL 命名提示它们可能存在平台、样本或批次差异，不能仅凭组织名称直接视为同分布。

### Breast 类别

- `Breast_Cancer`

只有一个主切片，不能在当前 7 个数据集中形成同组织训练/测试对照。

## 平台信息现状

当前本地文件能确认的只有部分平台线索：

- `HumanTonsil_10x`：名称明确包含 10x；
- `HumanTonsil_AddOns_S1`：名称明确包含 AddOns S1；
- `Breast_Cancer`：原始 registry 明确指向 CytAssist FFPE Gene + Protein Expression 文件；
- `A1_LN`、`D1_LN`、`A1_TNSL`、`D1_TNSL`：当前本地 registry 没有给出足够的正式平台/病人元数据，不能根据 H5AD 的字段名反推具体平台。

## 历史实验记录

历史的 3→4、4→3 和按组织整类留出结果已按当前实验记录清理，不作为项目现行结果保存。

## Downloads 压缩包中的补充证据

`C:\Users\xk2714\Downloads\ST2P (1).zip` 中的 `artifacts/metadata/scn_pro_final_matrix_20261023/` 保存了早期同组织实验目录。目录名称和结果 JSON 明确显示：

- `same_tissue_D1LN_to_A1LN`：用 `D1_LN` 训练/验证，测试 `A1_LN`；
- `same_tissue_A1LN_to_D1LN`：用 `A1_LN` 训练/验证，测试 `D1_LN`；
- `same_tissue_D1TNSL_to_A1TNSL`：用 `D1_TNSL` 训练/验证，测试 `A1_TNSL`；
- `same_tissue_A1TNSL_to_D1TNSL`：用 `A1_TNSL` 训练/验证，测试 `D1_TNSL`。

这说明原项目把 LN 和 TNSL 的同组织泛化定义为**成对切片之间的双向测试**，不是把 4 个 A1/D1 切片合成一个训练组、再把剩余 3 个切片作为另一组。

压缩包中还存在 `E1_Tonsil10x_to_AddOns` 和反向 `E1_reverse_AddOns_to_Tonsil10x`，它们被归档为 cross-slice 实验，不能自动等同于上面的 `same_tissue_*` 配对实验。

## 证据文件

- 原始注册表：`ST2P/data/raw/DGAT_final7_breast_raw_registry_20261021/registry_manifest.json`
- 冻结预处理清单：`ST2P/data/processed/DGAT_final7_breast_ptprc31_muon017_axis0_20261022/prepared_manifest.json`
- 预处理脚本：`ST2P/src/prepare_dgat_lodo_registry_muon_axis0.py`

官方仓库：

- `https://github.com/coderzhd/ST2P/blob/main/README.md`
- `https://github.com/coderzhd/ST2P/blob/main/src/build_final7_breast_raw_registry.py`
- `https://github.com/coderzhd/ST2P/blob/main/docs/REPRODUCIBILITY.md`

## 待补信息

需要从原始数据来源、论文补充材料或师兄处确认：

1. 每个数据集的准确组织和解剖部位；
2. 病人/样本编号；
3. 测序和空间检测平台；
4. 建库方式和批次；
5. A1/D1 与 `HumanTonsil_10x`/`AddOns_S1` 是否属于同一组织分布；
6. 哪些切片可以组成正式的同组织训练/测试组。
