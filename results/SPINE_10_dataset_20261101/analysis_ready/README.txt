SPINE 分析结果统一目录

每个实验一个文件夹；后续分析只读取本目录。
*_predictions.npy：预测蛋白表达矩阵
*_predictions.h5ad：带 spot barcode 和蛋白名称的分析文件
*_predictions.csv.gz：带行列名的压缩表格
*_truth.npy：真实蛋白表达
*_metrics.json / summary.json：评估指标

注意：Breast Cancer 跨组织方向已按组织级留出规则保存三个目标：`Breast_Cancer`、`PBC_PR_6837`、`PBC_PR_6835_5A`。
