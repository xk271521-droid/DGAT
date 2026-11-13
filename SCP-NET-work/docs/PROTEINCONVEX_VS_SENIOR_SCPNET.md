# ProteinConvex与师兄SCP-NET的实际代码对比

核对日期：2026-10-30。

比较基线是用户确认的师兄修改版 `ST2P/src/scn_pro.py`，
工作版本是 `experiments/protein_convex/scn_pro_protein_convex.py`。
两份当前源文件的SHA256均与各自保留的Breast运行脚本快照一致。
这里不把师兄文件等同于论文作者未核实的完整原仓库。

## 核心结论

师兄的多尺度主干和损失保留。累计新增两层机制：

1. CognateAnchor：把每个蛋白对应基因的三尺度RNA直接注入隐藏特征和预测残差。
2. ProteinConvex：在直接RNA残差通路中，用每蛋白一个alpha混合原动态尺度权重与静态可训练权重。

原动态多尺度隐藏特征融合不变；ProteinConvex的新混合发生在RNA锚定预测通路。

## 师兄已包含的机制

| 机制 | 师兄版本 | ProteinConvex |
|---|---|---|
| 共享RNA编码器和相对坐标编码 | 有 | 保留 |
| 中心/局部/区域三尺度表征 | 有 | 保留 |
| 蛋白特异邻居查询和注意力 | 有 | 保留 |
| 全局+低秩逐蛋白动态尺度路由 | 有 | 保留 |
| 每蛋白独立非线性预测头 | 有 | 保留 |
| RNA全切片上下文、受限校正、低秩上下文适配 | 有 | 保留 |
| 对应基因RNA显式锚定 | 无 | 新增 |
| RNA锚定通路的逐蛋白静态/动态混合 | 无 | 新增 |

## 新增1：对应基因RNA锚定

根据蛋白名查对应RNA基因。例如CD3E→CD3E，MS4A1→MS4A1，
PTPRC_1和PTPRC_2均映射到PTPRC。映射来自代码中的字典，没有引入外部训练数据。
当前保留checkpoint的mask显示30/31个蛋白启用，HLA_DRA未配置映射而关闭此分支；
这不等于断言数据中不存在HLA对应RNA。

对spot i和蛋白p，提取
`u_ip = [中心对应基因RNA, local4均值, regional4均值]`。
regional_neighbors=8在此代码表示总共8个邻居，前4个local，后4个regional。
u来自原输入路径，经过既有基因扰动和切片校正，不是额外的clean RNA流。
邻居均值直接对各自4个邻居求平均，没有根据attention token mask重新加权。

原context-adapted隐藏特征记作z_ip，原预测头记作f_p。
新增隐藏注入 `z'_ip = z_ip + A_p u_ip`，A_p是每蛋白的3→256线性映射，零初始化。
新增直接预测 `d_ip = w_p^T u_ip + b_p`。
最后 `yhat_ip = f_p(z'_ip) + g_p d_ip`，其中g_p=sigmoid(gate_logit_p)，
初始gate_logit=-2，对应g约0.119。g仅控制直接残差，A_p注入由A_p自身学习，
不能把g解释为同时控制整个锚定模块的门控。无映射时mask关闭两条新增路径。
隐藏注入从0开始，但直接支路初始非零，不能宣称初始预测严格等价于无锚定模型。

锚定接到final及center/local/regional辅助预测分支，原辅助监督仍在。

## 新增2：ProteinConvex在锚定通路混合尺度

师兄已有每spot、每蛋白动态尺度向量a_ip=softmax(router(...))。
CognateAnchor中的w_p是每蛋白一组、跨spots共享的可训练静态权重，初始[0.5,0.3,0.2]。
ProteinConvex新增alpha_p=sigmoid(alpha_logit_p)，每蛋白一个，跨spots共享，初始约0.119。

`v_ip = alpha_p * a_ip + (1-alpha_p) * w_p`

`d_ip = v_ip^T u_ip + b_p`

其余隐藏注入和最终预测公式沿用CognateAnchor。
alpha_p决定直接锚定预测多大程度使用spot动态尺度偏好，多少使用该蛋白的静态尺度偏好。
动态a_ip随spot变化；alpha_p自身不是由每个spot额外计算的动态门。
只有final分支传入a_ip；三个辅助分支仍使用静态w_p，不执行这次混合。

## “Convex”的准确边界

alpha_p在0到1，因此v_ip是两套权重向量之间的有界插值。
但是原版direct_weight是无约束Parameter，没有softmax，不保证三项非负或和为1。
初值[0.5,0.3,0.2]的和是1，不代表训练后仍为1。
保留Breast checkpoint中静态权重行和最大约0.7118，确实不满足概率单纯形约束。
因此不能把当前实现描述成“严格的三尺度概率凸组合”，也不能把混合只理解为概率分配：
它还可能改变直接RNA残差的总体幅度。
后来单独尝试的静态softmax版本不是当前保留的ProteinConvex版本。

## 新增参数

| 新增参数 | 数量 |
|---|---:|
| A_p：31×3×256 | 23,808 |
| 静态尺度权重31×3 | 93 |
| 直接偏置31 | 31 |
| 直接残差gate31 | 31 |
| ProteinConvex的alpha31 | 31 |
| 合计 | 23,994 |

师兄4,706,592 → CognateAnchor4,730,555 → ProteinConvex4,730,586。
累计参数增加约0.51%；从CognateAnchor到ProteinConvex只新增31个alpha参数。

## 训练设置与实验流程

两者保留的配置中，7切片注册表、9764基因/31蛋白、Muon CLR axis=0、
6源1目标、seed42、20epochs、batch64、0.1验证比例、lr0.0005、weight_decay0.0001、
邻居local4/regional4、gene_dropout0.5、neighbor_mask0.25等相同。
损失均为 `MSE_final + 0.25*MSE_aux_mean + 0.05*PCC_loss`，没有新增锚定专用损失。

实验流程另有收紧：师兄文件每epoch计算目标诊断，并另存test-oracle checkpoint；
正式source-selected checkpoint依然由源验证选出。ProteinConvex去掉逐epoch目标诊断和oracle分支，
目标ADT数值仅在训练及选点后读取，面板名称元数据仍按原流程读取。
不能把这个流程修改描述成“师兄正式结果是用目标选点”，以下对照均取source-selected。

## 保留的Breast结果

| 方法 | PCC-Pro ↑ | PCC-Spot ↑ | RMSE-Global ↓ |
|---|---:|---:|---:|
| 师兄scn_pro.py | 0.5952222576 | 0.5926463523 | 0.3130912781 |
| CognateAnchor | 0.5954097496 | 0.6286801792 | 0.3018039167 |
| ProteinConvex | 0.6020346726 | 0.6332672906 | 0.2995744348 |

这说明保留结果中三项均有改善，不构成单seed下每个新增机制独立因果贡献的证明。
原始证据分别在 `runs/baseline_Breast_Cancer/results.json`、
`runs/strict_cognate_anchor_Breast_Cancer/results.json`、
`runs/protein_convex_Breast_Cancer/results.json`。

## 主要代码位置

- 师兄原有尺度融合及预测：`ST2P/src/scn_pro.py:204`。
- 新锚定模块：`experiments/protein_convex/scn_pro_protein_convex.py:91`。
- alpha及静态/动态混合：同文件103、117行。
- u提取与解码接入：同文件252行。
- 对应基因字典：同文件483行。
- 原损失与源验证选点：同文件538、544行。
- 目标最终读取：同文件556行。
