# 快速开始 - Mast Cell 数据训练

## Windows 用户

### 1. 准备数据 (在 R 中)

打开 R 或 RStudio，运行：
```r
setwd("E:/Development/Squidiff")
source("scripts/1.data_output.r")
```

### 2. 验证数据维度
```bash
python E:\Development\Squidiff\scripts\check_shape.py --data_path "E:\Development\Squidiff\data\mast_cells_200genes.h5ad"
# Data shape: (402, 159): (402, 159)
```

### 3. 填写数据维度开始炼丹
```bash
# train mast cells
## genesize和outputdim填写Data shape: (402, 159)的159
python train_squidiff.py --logger_path "./results" --data_path "E:\Development\Squidiff\data\mast_cells_200genes.h5ad" --resume_checkpoint "./models" --gene_size 159 --output_dim 159
```
