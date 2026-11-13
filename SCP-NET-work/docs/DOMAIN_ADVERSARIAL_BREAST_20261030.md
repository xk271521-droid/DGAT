# ProteinConvex域对抗表示：失败记录

日期：2026-10-30。候选从ProteinConvex派生，在RNA-only slice context上加入6类源切片分类器和gradient reversal。
domain loss权重0.02、GRL系数1.0，其他参数、数据、mask、20 epoch和source-only选点不变。

结果：源验证最佳epoch18，PCC-Pro=0.752573；Breast测试PCC-Pro=0.556399、PCC-Spot=0.544199、RMSE=0.334138。
相比ProteinConvex 0.602035/0.633267/0.299574，变化为-0.045636/-0.089068/+0.034564，三项均明显退化。
域分类准确率约0.34–0.43，接近6类随机基线0.167但仍能识别部分源差异，说明约束实际生效但破坏了当前蛋白预测表示。

**决策：拒绝。** 不扩展七折，不保留候选代码、checkpoint、预测或运行日志；只保留本失败记录。ProteinConvex正式基线不变。

清理日期：2026-10-30。已删除：

- `experiments/protein_convex_domain_adversarial_20261030/`
- `runs/protein_convex_domain_adversarial_Breast_Cancer_20261030/`
- `runs/protein_convex_domain_adversarial_Breast_Cancer_20261030.log`

该方法不是DGAT或DPAS-Graph官方代码复现，也没有修改原作者仓库、冻结数据或正式汇总表。
