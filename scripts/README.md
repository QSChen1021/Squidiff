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

### 4. 运行推理
```bash
python scripts\3.sample_mast_cells.py ^
  --model_path logs\mast_cells_<时间戳>\model.pt ^
  --data_path data\mast_cells_200genes.h5ad ^
  --output_dir outputs
```

---

## Linux/Mac 用户

### 1. 准备数据 (在 R 中)
```r
setwd("E:/Development/Squidiff")
source("scripts/1.data_output.r")
```

### 2. 格式转换 (Python)
```bash
cd E:/Development/Squidiff
python scripts/5.convert_h5seurat_to_h5ad.py data/mast_cells_200genes.h5seurat
```

### 3. 验证数据
```bash
python scripts/4.validate_data.py --data_path data/mast_cells_200genes.h5ad
```

### 4. 训练模型
```bash
bash scripts/2.train_mast_cells.sh
```

### 5. 运行推理
```bash
python scripts/3.sample_mast_cells.py \
  --model_path logs/mast_cells_<时间戳>/model.pt \
  --data_path data/mast_cells_200genes.h5ad \
  --output_dir outputs
```

---

## 输出文件

| 文件 | 位置 | 说明 |
|------|------|------|
| R 输出 | `data/mast_cells_200genes.h5seurat` | R 脚本输出 (中间格式) |
| 数据文件 | `data/mast_cells_200genes.h5ad` | Python 转换后 (最终格式) |
| 模型文件 | `logs/mast_cells_*/model.pt` | 训练好的模型 |
| 预测结果 | `outputs/mast_cells_predicted.h5ad` | 模型预测结果 |

---

## 故障排除

### 问题 1: 文件格式不匹配

如果验证脚本提示 `.h5seurat` 文件：
```bash
python scripts/5.convert_h5seurat_to_h5ad.py data/mast_cells_200genes.h5seurat
```

### 问题 2: 找不到 data 目录
```bash
mkdir data
```

### 问题 3: 找不到 logs 目录
```bash
mkdir logs
```

### 问题 3: CUDA out of memory
编辑 `scripts/2.train_mast_cells.bat`，减小 BATCH_SIZE:
```batch
set BATCH_SIZE=32
```

### 问题 4: 基因数不匹配错误
检查数据实际基因数，然后更新脚本中的 GENE_SIZE 参数。
