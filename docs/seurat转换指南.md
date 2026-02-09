# Seurat 数据转换指南

> 将 Seurat (R) 对象转换为 Squidiff 所需的 h5ad 格式

## 目录

- [概述](#概述)
- [环境准备](#环境准备)
- [转换流程](#转换流程)
- [数据格式要求](#数据格式要求)
- [完整示例](#完整示例)
- [常见问题](#常见问题)

---

## 概述

Squidiff 使用 Python 的 AnnData 格式（.h5ad 文件）存储单细胞数据。如果您使用的是 Seurat (R) 进行数据分析，需要先将数据转换为 h5ad 格式。

本指南提供两种转换方法：

1. **使用 `SeuratDisk`**（推荐）- 专门的 R/Python 数据转换工具
2. **使用 `scRNA-seq`** - 传统的单细胞数据转换流程

---

## 环境准备

### 方法一：使用 SeuratDisk（推荐）

#### R 环境（导出 Seurat 对象）

```r
# 安装 SeuratDisk
if (!require("remotes")) install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")

# 加载库
library(Seurat)
library(SeuratDisk)
```

#### Python 环境（导入并转换为 h5ad）

```bash
# 安装依赖
uv pip install anndata scipy scikit-learn
```

### 方法二：使用 loom/py文件

```r
# 安装依赖
install.packages("Seurat")
install.packages("reticulate")
```

---

## 转换流程

### 方法一：SeuratDisk 转换（推荐）

#### 步骤 1：在 R 中保存 Seurat 对象为 .h5seurat 格式

```r
# 加载您的 Seurat 对象
# 假设您的对象名为 seurat_obj
library(Seurat)
library(SeuratDisk)

# 检查对象结构
print(seurat_obj)
# 检查元数据
head(seurat_obj@meta.data)

# 保存为 .h5seurat 格式
Save(seurat_obj, filename = "my_data.h5seurat")
```

#### 步骤 2：在 Python 中转换为 .h5ad 格式

```python
from seurat_disk import SeuratDisk

# 转换为 h5ad
SeuratDisk.convert(
    "my_data.h5seurat",
    dest = "h5ad",
    output = "my_data.h5ad"
)

# 验证转换结果
import scanpy as sc
adata = sc.read_h5ad("my_data.h5ad")
print(adata)
```

### 方法二：通过 Loom 格式转换

#### 步骤 1：在 R 中导出为 Loom 格式

```r
library(Seurat)

# 确保 Seurat 对象已经标准化
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# 导出为 loom 格式
library(Seurat)
library(loomR)

# 创建 loom 文件
lfile <- create("my_data.loom")
# 添加矩阵
lfile[["matrix"]] <- as.matrix(GetAssayData(seurat_obj, slot = "data"))
# 添加元数据
lfile$col.attrs$obs <- seurat_obj@meta.data
# 添加基因信息
lfile$row.attrs$var <- data.frame(
  gene_names = rownames(seurat_obj),
  row.names = rownames(seurat_obj)
)
# 保存并关闭
close(lfile)
```

#### 步骤 2：在 Python 中读取 Loom 并保存为 h5ad

```python
import scanpy as sc
import anndata

# 读取 loom 文件
adata = sc.read_loom("my_data.loom")

# 保存为 h5ad
adata.write_h5ad("my_data.h5ad")
```

### 方法三：直接导出 CSV/H5（适用于小数据集）

#### 步骤 1：在 R 中导出

```r
library(Seurat)

# 导出表达矩阵
write.csv(as.matrix(GetAssayData(seurat_obj, slot = "data")),
          file = "expression_matrix.csv")

# 导出元数据
write.csv(seurat_obj@meta.data, file = "metadata.csv")
```

#### 步骤 2：在 Python 中导入并创建 AnnData 对象

```python
import pandas as pd
import scanpy as sc
from anndata import AnnData
import numpy as np

# 读取数据
expr_df = pd.read_csv("expression_matrix.csv", index_col=0)
meta_df = pd.read_csv("metadata.csv", index_col=0)

# 创建 AnnData 对象
adata = AnnData(
    X=expr_df.values,
    obs=meta_df,
    var=pd.DataFrame(index=expr_df.columns)
)

# 保存为 h5ad
adata.write_h5ad("my_data.h5ad")
```

---

## 数据格式要求

Squidiff 需要以下数据结构：

### 1. 基础数据结构

```python
import scanpy as sc

# 必需字段
adata.X                    # 基因表达矩阵 (cells x genes)
adata.obs                  # 细胞元数据 (DataFrame)
adata.var                  # 基因信息 (DataFrame)

# 示例数据结构
# adata.obs 包含以下列：
# - Group: 细胞分组/条件（必需）
# - cell_type: 细胞类型（可选）
# - batch: 批次信息（可选）
```

### 2. 药物扰动数据（使用药物结构时）

如果使用 `--use_drug_structure True`，需要在元数据中添加：

```python
# 在 adata.obs 中添加：
adata.obs['SMILES'] = [
    'CC(=O)OC1=CC=CC=C1C(=O)O',  # 阿司匹林示例
    'CC(=O)NC1=CC=C(C=C1)O',     # 对乙酰氨基酚示例
    # ... 每个细胞对应的药物 SMILES
]

adata.obs['dose'] = [
    1.0,   # 药物剂量（μM 或其他单位）
    0.5,
    # ... 每个细胞对应的药物剂量
]
```

### 3. 对照组数据

使用药物结构时，需要单独的对照组数据：

```python
# 对照组 AnnData 对象，不包含药物信息
control_adata = adata[adata.obs['condition'] == 'control']
control_adata.write_h5ad("control_data.h5ad")
```

---

## 完整示例

### 示例 1：药物扰动数据转换

#### R 端代码

```r
library(Seurat)
library(SeuratDisk)

# 假设您有一个药物处理的 Seurat 对象
# 确保元数据包含以下列：
# - condition: 处理条件（如 'drug_A', 'control'）
# - drug_name: 药物名称
# - dose: 药物剂量

# 添加 SMILES 列（如果还没有）
smiles_map <- c(
  "drug_A" = "CC(=O)OC1=CC=CC=C1C(=O)O",  # 阿司匹林
  "drug_B" = "CC(=O)NC1=CC=C(C=C1)O"      # 对乙酰氨基酚
)

seurat_obj$SMILES <- smiles_map[seurat_obj$drug_name]

# 检查数据
head(seurat_obj@meta.data)

# 保存为 h5seurat
Save(seurat_obj, filename = "drug_perturbation.h5seurat")
```

#### Python 端代码

```python
from seurat_disk import SeuratDisk
import scanpy as sc
import pandas as pd

# 转换为 h5ad
SeuratDisk.convert(
    "drug_perturbation.h5seurat",
    dest="h5ad",
    output="drug_perturbation.h5ad"
)

# 读取并验证
adata = sc.read_h5ad("drug_perturbation.h5ad")
print(adata)
print(adata.obs.columns)

# 检查必需的列
assert 'SMILES' in adata.obs.columns, "缺少 SMILES 列"
assert 'dose' in adata.obs.columns, "缺少 dose 列"
assert 'Group' in adata.obs.columns, "缺少 Group 列"

# 分离对照组
control_mask = adata.obs['condition'] == 'control'
control_adata = adata[control_mask].copy()
treated_adata = adata[~control_mask].copy()

# 保存
control_adata.write_h5ad("control_data.h5ad")
treated_adata.write_h5ad("treated_data.h5ad")

print("转换完成！")
print(f"对照组细胞数: {control_adata.n_obs}")
print(f"处理组细胞数: {treated_adata.n_obs}")
```

### 示例 2：细胞分化轨迹数据

#### R 端代码

```r
library(Seurat)
library(SeuratDisk)

# 假设您有时间序列的分化数据
# 元数据包含：
# - time_point: 时间点（如 'day0', 'day1', 'day5'）
# - cell_type: 细胞类型

# 创建 Group 列（Squidiff 需要此列）
seurat_obj$Group <- paste0(
  seurat_obj$time_point,
  "_",
  seurat_obj$cell_type
)

# 保存
Save(seurat_obj, filename = "differentiation.h5seurat")
```

#### Python 端代码

```python
from seurat_disk import SeuratDisk
import scanpy as sc

# 转换
SeuratDisk.convert(
    "differentiation.h5seurat",
    dest="h5ad",
    output="differentiation.h5ad"
)

# 读取和验证
adata = sc.read_h5ad("differentiation.h5ad")
print(adata.obs['Group'].value_counts())

# 数据预处理（可选）
import scvelo as scv
# 进行质量控制和标准化
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.raw = adata.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 保存处理后的数据
adata.write_h5ad("differentiation_processed.h5ad")
```

---

## 常见问题

### Q1: 转换后数据矩阵为空或全为零

**原因**: 可能是使用了错误的 slot（如 `counts` 而非 `data`）

**解决**: 在 R 中确保导出正确的数据

```r
# 使用标准化后的数据
GetAssayData(seurat_obj, slot = "data")  # 推荐

# 或使用原始计数
GetAssayData(seurat_obj, slot = "counts")
```

### Q2: 元数据丢失或格式错误

**原因**: SeuratDisk 版本不兼容或元数据格式问题

**解决**: 检查并转换元数据格式

```r
# 确保 meta.data 是 data.frame
seurat_obj@meta.data <- as.data.frame(seurat_obj@meta.data)

# 检查列名
print(colnames(seurat_obj@meta.data))
```

### Q3: 基因名称不匹配

**原因**: R 和 Python 对基因名的处理可能不同

**解决**: 统一基因名格式

```python
# 在 Python 中处理基因名
adata.var_names = adata.var_names.str.upper()  # 转为大写
# 或
adata.var_names = adata.var_names.str.strip()  # 去除空格
```

### Q4: 内存不足

**原因**: 数据集过大

**解决**: 分批处理或使用稀疏矩阵

```r
# 在 R 中创建稀疏矩阵
library(Matrix)
seurat_obj@assays$RNA@data <- Matrix(
  seurat_obj@assays$RNA@data,
  sparse = TRUE
)
```

### Q5: SMILES 格式错误

**原因**: SMILES 字符串格式不正确

**解决**: 验证 SMILES 格式

```python
from rdkit import Chem

# 验证 SMILES
def validate_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False

# 检查所有 SMILES
adata.obs['valid_smiles'] = adata.obs['SMILES'].apply(validate_smiles)
print(adata.obs['valid_smiles'].value_counts())

# 移除无效的 SMILES
adata = adata[adata.obs['valid_smiles']].copy()
```

---

## 附录：快速参考

### R 端常用命令

```r
# 查看 Seurat 对象
print(seurat_obj)
head(seurat_obj@meta.data)
colnames(seurat_obj@meta.data)

# 添加新列到元数据
seurat_obj$new_column <- values

# 筛选细胞
seurat_obj_subset <- subset(seurat_obj, condition == "treated")

# 保存
Save(seurat_obj, filename = "data.h5seurat")
```

### Python 端常用命令

```python
import scanpy as sc

# 读取数据
adata = sc.read_h5ad("data.h5ad")

# 查看数据
print(adata)
adata.obs.head()
adata.var.head()

# 筛选细胞
adata_subset = adata[adata.obs['condition'] == 'treated'].copy()

# 保存
adata.write_h5ad("data.h5ad")
```

---

## 相关资源

- **Seurat 文档**: https://satijalab.org/seurat/
- **SeuratDisk 文档**: https://mojaveazure.github.io/seurat-disk/
- **AnnData 文档**: https://anndata.readthedocs.io/
- **Scanpy 文档**: https://scanpy.readthedocs.io/

---

*文档最后更新：2025年*
