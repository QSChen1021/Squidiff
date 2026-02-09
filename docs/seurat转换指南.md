# Seurat 数据转换指南

> 将 Seurat (R) 对象转换为 Squidiff 所需的 h5ad 格式

## 目录

- [概述](#概述)
- [数据要求详解](#数据要求详解)
- [环境准备](#环境准备)
- [转换流程](#转换流程)
- [完整示例](#完整示例)
- [常见问题](#常见问题)

---

## 概述

Squidiff 使用 Python 的 AnnData 格式（.h5ad 文件）存储单细胞数据。如果您使用的是 Seurat (R) 进行数据分析，需要先将数据转换为 h5ad 格式。

本指南提供：
1. **详细的数据要求说明** - 包括细胞数、基因数、metadata 列等
2. **完整的转换流程** - 从 Seurat 到 h5ad 的多种方法
3. **验证检查清单** - 确保数据符合 Squidiff 要求

---

## 数据要求详解

### 核心数据结构

Squidiff 基于 AnnData 对象存储数据，核心结构如下：

| 组件 | 说明 | 数据类型 | 必需 |
|------|------|----------|------|
| `adata.X` | 基因表达矩阵 | `numpy.ndarray` 或 `scipy.sparse.spmatrix` | ✅ 是 |
| `adata.obs` | 细胞元数据 | `pandas.DataFrame` | ✅ 是 |
| `adata.var` | 基因信息 | `pandas.DataFrame` | ✅ 是 |
| `adata.obs['Group']` | 细胞分组标签 | `str`/`int` | ✅ 是 |

### 维度要求

#### 1. 基因数 (Gene Size)

基因数必须与训练时的 `--gene_size` 参数**完全一致**。

```python
# 示例：如果训练时设置 gene_size=500
adata.n_vars  # 必须等于 500
```

**常见基因数范围**：
- **小规模**: 100-500 基因（快速测试）
- **中规模**: 500-2000 基因（常用）
- **大规模**: 2000-10000 基因（完整转录组）

**重要**: `gene_size` 既影响输入维度，也影响模型结构参数。

#### 2. 细胞数 (Cell Count)

| 细胞数范围 | 推荐场景 | batch_size 建议 |
|-----------|---------|----------------|
| 1,000 - 10,000 | 小型数据集/快速测试 | 32-64 |
| 10,000 - 100,000 | 标准数据集 | 64-128 |
| >100,000 | 大规模数据集 | 128-256 |

**最低要求**: 建议至少 **500 个细胞**以确保模型训练稳定性。

#### 3. 表达矩阵格式

```python
# Squidiff 支持两种格式：

# 格式 1: 密集矩阵 (numpy.ndarray)
type(adata.X) == np.ndarray  # True

# 格式 2: 稀疏矩阵 (scipy.sparse)
from scipy.sparse import issparse
issparse(adata.X)  # True
```

**推荐**: 对于超过 10,000 基因的数据集，使用**稀疏矩阵**以节省内存。

### Metadata (元数据) 要求

#### 必需列

| 列名 | 数据类型 | 说明 | 使用场景 |
|------|---------|------|----------|
| `Group` | `str` 或 `int` | 细胞分组/条件标签 | **所有场景必需** |

#### 使用药物结构时的必需列

当 `--use_drug_structure True` 时，需要额外添加：

| 列名 | 数据类型 | 说明 | 示例值 |
|------|---------|------|--------|
| `SMILES` | `str` | 药物 SMILES 字符串 | `"CC(=O)OC1=CC=CC=C1C(=O)O"` |
| `dose` | `float` | 药物剂量 | `1.0` (μM) |

**SMILES 格式要求**：
- 必须是有效的 SMILES 字符串
- 药物组合用 `+` 连接（如 `comb_num > 1`）
- 示例：`"CC(=O)OC1=CC=CC=C1C(=O)O+CC(=O)NC1=CC=C(C=C1)O"`

#### 可选列（推荐添加）

| 列名 | 说明 | 用途 |
|------|------|------|
| `cell_type` | 细胞类型 | 结果分析和可视化 |
| `condition` | 实验条件 | 对照组/处理组筛选 |
| `batch` | 批次信息 | 批次效应校正 |
| `time_point` | 时间点 | 时间序列分析 |

### 对照组数据要求

使用药物结构时 (`--use_drug_structure True`)：

1. **必需单独的对照组文件**
2. 对照组细胞的 `gene_size` 必须与处理组**完全一致**
3. 对照组**不需要** `SMILES` 和 `dose` 列

```python
# 对照组数据结构示例
control_adata.X.shape      # (N_cells, gene_size)
'SMILES' not in control_adata.obs.columns  # True - 对照组不需要
```

### 数据预处理建议

#### 在 R 中预处理（推荐）

```r
# 标准单细胞预处理流程
library(Seurat)

# 1. 质量控制
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)
seurat_obj <- subset(seurat_obj, subset = percent.mt < 20)

# 2. 标准化
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 3. 如果使用特定基因集
# seurat_obj <- seurat_obj[VariableFeatures(seurat_obj), ]
```

#### 在 Python 中预处理

```python
import scanpy as sc

# 质量控制
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 线性回归校正线粒体基因比例
adata.obs['pct_counts_mt'] = np.sum(
    adata[:, adata.var_names.str.startswith('MT-')].X, axis=1
) / np.sum(adata.X, axis=1) * 100

# 标准化
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

---

## 环境准备

### 方法一：使用 SeuratDisk（推荐 - R 端直接导出 h5ad）

**注意**: SeuratDisk 是一个 **R 包**，需要在 R 中使用。

#### R 环境（导出 Seurat 对象为 h5ad）

```r
# 安装 SeuratDisk
if (!require("remotes")) install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")

# 加载库
library(Seurat)
library(SeuratDisk)
```

#### Python 环境（直接读取 h5ad）

```bash
# 安装依赖
uv pip install anndata scipy scikit-learn scanpy h5py
```

**重要**: SeuratDisk 可以在 R 中直接将 Seurat 对象保存为 `.h5ad` 格式，Python 端只需用 `scanpy` 读取即可。

### 方法二：使用 loom/csv 中间格式

```r
# 安装依赖
install.packages("Seurat")
install.packages("reticulate")
```

---

## 转换流程

### 方法一：SeuratDisk 直接导出 h5ad（推荐）

**说明**: SeuratDisk 支持**直接导出为 h5ad 格式**，无需 Python 端转换。

#### 步骤 1：在 R 中直接保存为 .h5ad 格式

```r
library(Seurat)
library(SeuratDisk)

# 检查对象结构
print(seurat_obj)
print(paste("细胞数:", ncol(seurat_obj)))
print(paste("基因数:", nrow(seurat_obj)))

# 检查元数据
head(seurat_obj@meta.data)

# 确保 Group 列存在
if (!"Group" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$Group <- seurat_obj$seurat_clusters  # 或其他分组信息
}

# 保存为 .h5ad 格式（SeuratDisk 支持）
SaveH5Seurat(seurat_obj, filename = "my_data.h5ad")
```

#### 步骤 2：在 Python 中直接读取

```python
import scanpy as sc

# 直接读取 h5ad 文件
adata = sc.read_h5ad("my_data.h5ad")

# 验证
print(f"细胞数: {adata.n_obs}")
print(f"基因数: {adata.n_vars}")
print(f"元数据列: {adata.obs.columns.tolist()}")
```

**注意**: 如果 `SaveH5Seurat` 函数不可用，可以使用以下替代方法：

```r
# 方法 A: 保存为 .h5seurat，然后在 R 中转换
Save(seurat_obj, filename = "my_data.h5seurat")
Convert("my_data.h5seurat", dest = "h5ad", overwrite = TRUE)

# 方法 B: 使用 reticulate 调用 Python 的 anndata
library(reticulate)
anndata <- import("anndata")
adata <- anndata$AnnData(
  X = as.matrix(GetAssayData(seurat_obj, slot = "data")),
  obs = seurat_obj@meta.data,
  var = data.frame(row.names = rownames(seurat_obj))
)
adata$write_h5ad("my_data.h5ad")
```

### 方法二：通过 Loom 格式转换

#### 步骤 1：在 R 中导出为 Loom 格式

```r
library(Seurat)
library(loomR)

# 确保 Seurat 对象已经标准化
seurat_obj <- NormalizeData(seurat_obj)

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
from anndata import AnnData

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

## 完整示例

### 示例 1：药物扰动数据转换

#### R 端代码

```r
library(Seurat)
library(SeuratDisk)

# ========== 数据检查 ==========
print(paste("细胞数:", ncol(seurat_obj)))
print(paste("基因数:", nrow(seurat_obj)))
print(colnames(seurat_obj@meta.data))

# ========== 准备元数据 ==========
# 确保 Group 列存在（Squidiff 必需）
if (!"Group" %in% colnames(seurat_obj@meta.data)) {
  # 创建 Group 列：条件_细胞类型
  seurat_obj$Group <- paste0(
    seurat_obj$condition,
    "_",
    seurat_obj$cell_type
  )
}

# 添加 SMILES 列（使用药物结构时必需）
smiles_map <- c(
  "drug_A" = "CC(=O)OC1=CC=CC=C1C(=O)O",  # 阿司匹林
  "drug_B" = "CC(=O)NC1=CC=C(C=C1)O",     # 对乙酰氨基酚
  "drug_C" = "CC(C)CC1=CC=C(C=C1)C(C)C=C2"  # 布洛芬
)

seurat_obj$SMILES <- smiles_map[seurat_obj$drug_name]

# 检查 dose 列
if (!"dose" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$dose <- 1.0  # 默认剂量
}

# 验证元数据
required_cols <- c("Group", "SMILES", "dose")
missing_cols <- required_cols[!required_cols %in% colnames(seurat_obj@meta.data)]
if (length(missing_cols) > 0) {
  stop(paste("缺少必需列:", paste(missing_cols, collapse=", ")))
}

# ========== 基因选择（可选） ==========
# 如果需要限制基因数，选择高变基因
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
variable_genes <- VariableFeatures(seurat_obj)
seurat_obj_subset <- seurat_obj[variable_genes, ]

print(paste("选择的基因数:", nrow(seurat_obj_subset)))

# ========== 保存 ==========
Save(seurat_obj_subset, filename = "drug_perturbation.h5seurat")
```

#### Python 端代码

```python
import scanpy as sc
import pandas as pd
import numpy as np

# ========== 直接读取 h5ad（已在 R 端转换） ==========
adata = sc.read_h5ad("drug_perturbation.h5ad")

print(f"细胞数: {adata.n_obs}")
print(f"基因数: {adata.n_vars}")
print(f"数据类型: {type(adata.X)}")
print(f"元数据列: {adata.obs.columns.tolist()}")

# ========== 检查必需列 ==========
required_cols = ['Group', 'SMILES', 'dose']
for col in required_cols:
    if col not in adata.obs.columns:
        raise ValueError(f"缺少必需列: {col}")
    print(f"✓ {col}: {adata.obs[col].nunique()} 个唯一值")

# ========== 检查 SMILES 格式 ==========
from rdkit import Chem

def validate_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False

adata.obs['valid_smiles'] = adata.obs['SMILES'].apply(validate_smiles)
print(f"\n有效 SMILES: {adata.obs['valid_smiles'].sum()}/{len(adata)}")

if adata.obs['valid_smiles'].sum() < len(adata):
    print("⚠️ 警告: 存在无效的 SMILES 字符串")
    invalid = adata[~adata.obs['valid_smiles']]
    print(invalid.obs[['SMILES', 'drug_name']])

# ========== 分离对照组和处理组 ==========
if 'condition' in adata.obs.columns:
    control_mask = adata.obs['condition'] == 'control'
    control_adata = adata[control_mask].copy()
    treated_adata = adata[~control_mask].copy()

    print(f"\n对照组细胞数: {control_adata.n_obs}")
    print(f"处理组细胞数: {treated_adata.n_obs}")

    # 保存
    control_adata.write_h5ad("control_data.h5ad")
    treated_adata.write_h5ad("treated_data.h5ad")
else:
    # 如果没有 condition 列，保存全部数据
    adata.write_h5ad("drug_data.h5ad")
```

### 示例 2：细胞分化轨迹数据

#### R 端代码

```r
library(Seurat)
library(SeuratDisk)

# ========== 检查数据 ==========
print(paste("细胞数:", ncol(seurat_obj)))
print(paste("基因数:", nrow(seurat_obj)))

# ========== 创建 Group 列 ==========
# Squidiff 需要 Group 列来区分不同条件
seurat_obj$Group <- paste0(
  seurat_obj$time_point,
  "_",
  seurat_obj$cell_type
)

# 查看分组情况
print(table(seurat_obj$Group))

# ========== 保存 ==========
Save(seurat_obj, filename = "differentiation.h5seurat")
```

#### Python 端代码

```python
import scanpy as sc

# ========== 直接读取 h5ad（已在 R 端转换） ==========
adata = sc.read_h5ad("differentiation.h5ad")

print(f"细胞数: {adata.n_obs}")
print(f"基因数: {adata.n_vars}")
print("\nGroup 分布:")
print(adata.obs['Group'].value_counts())

# ========== 数据预处理（可选） ==========
# 质量控制
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 标准化
adata.raw = adata.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ========== 保存 ==========
adata.write_h5ad("differentiation_processed.h5ad")

print("\n✓ 转换完成！")
print(f"最终细胞数: {adata.n_obs}")
print(f"最终基因数: {adata.n_vars}")
```

### 示例 3：数据验证函数

在 Python 中创建一个完整的验证函数：

```python
import scanpy as sc
from rdkit import Chem
import numpy as np

def validate_squidiff_data(adata, use_drug_structure=False):
    """
    验证 AnnData 对象是否符合 Squidiff 要求

    参数:
        adata: AnnData 对象
        use_drug_structure: 是否使用药物结构

    返回:
        is_valid: bool
        issues: list of str
    """
    issues = []

    # 检查基本结构
    if adata.X is None:
        issues.append("❌ 缺少表达矩阵 (adata.X)")
        return False, issues

    # 检查维度
    n_cells, n_genes = adata.n_obs, adata.n_vars
    print(f"细胞数: {n_cells}")
    print(f"基因数: {n_genes}")

    if n_cells < 500:
        issues.append(f"⚠️ 细胞数较少 ({n_cells})，建议至少 500 个")

    if n_genes < 50:
        issues.append(f"⚠️ 基因数过少 ({n_genes})，建议至少 100 个")

    # 检查必需列
    if 'Group' not in adata.obs.columns:
        issues.append("❌ 缺少必需列: Group")
    else:
        print(f"✓ Group 列存在，{adata.obs['Group'].nunique()} 个唯一值")

    # 检查药物结构相关列
    if use_drug_structure:
        if 'SMILES' not in adata.obs.columns:
            issues.append("❌ 使用药物结构时缺少必需列: SMILES")
        else:
            # 验证 SMILES 格式
            valid_count = 0
            for smiles in adata.obs['SMILES']:
                if Chem.MolFromSmiles(smiles) is not None:
                    valid_count += 1
            if valid_count < len(adata):
                issues.append(f"⚠️ 有 {len(adata) - valid_count} 个无效的 SMILES")
            else:
                print(f"✓ SMILES 格式正确")

        if 'dose' not in adata.obs.columns:
            issues.append("❌ 使用药物结构时缺少必需列: dose")
        else:
            print(f"✓ dose 列存在，范围: {adata.obs['dose'].min():.2f} - {adata.obs['dose'].max():.2f}")

    # 检查数据类型
    if isinstance(adata.X, np.ndarray):
        print(f"✓ 密集矩阵格式")
    else:
        from scipy.sparse import issparse
        if issparse(adata.X):
            print(f"✓ 稀疏矩阵格式")
        else:
            issues.append(f"⚠️ 未知的数据类型: {type(adata.X)}")

    # 检查缺失值
    if np.isnan(adata.X.data if hasattr(adata.X, 'data') else adata.X).any():
        issues.append("⚠️ 存在缺失值 (NaN)")

    # 总结
    if len(issues) == 0:
        print("\n✓ 数据验证通过！")
        return True, []
    else:
        print("\n发现问题:")
        for issue in issues:
            print(issue)
        return False, issues

# 使用示例
adata = sc.read_h5ad("my_data.h5ad")
is_valid, issues = validate_squidiff_data(adata, use_drug_structure=True)
```

---

## 常见问题

### Q1: 转换后数据矩阵为空或全为零

**原因**: 可能是使用了错误的 slot（如 `counts` 而非 `data`）

**解决**:
```r
# 使用标准化后的数据
GetAssayData(seurat_obj, slot = "data")  # 推荐

# 或使用原始计数
GetAssayData(seurat_obj, slot = "counts")
```

### Q2: 元数据丢失或格式错误

**原因**: SeuratDisk 版本不兼容或元数据格式问题

**解决**:
```r
# 确保 meta.data 是 data.frame
seurat_obj@meta.data <- as.data.frame(seurat_obj@meta.data)

# 检查列名
print(colnames(seurat_obj@meta.data))
```

### Q3: 基因数不匹配错误

**错误信息**: `RuntimeError: mat1 and mat2 shapes cannot be multiplied`

**原因**: `gene_size` 参数与实际数据基因数不一致

**解决**:
```python
# 检查实际基因数
adata.n_vars

# 训练时使用正确的 gene_size
python train_squidiff.py --gene_size {adata.n_vars} ...
```

### Q4: 基因名称不匹配

**原因**: R 和 Python 对基因名的处理可能不同

**解决**:
```python
# 统一基因名格式
adata.var_names = adata.var_names.str.upper()  # 转为大写
# 或
adata.var_names = adata.var_names.str.strip()  # 去除空格
```

### Q5: 内存不足

**原因**: 数据集过大

**解决**:
```r
# 在 R 中创建稀疏矩阵
library(Matrix)
seurat_obj@assays$RNA@data <- Matrix(
  seurat_obj@assays$RNA@data,
  sparse = TRUE
)
```

### Q6: SMILES 格式错误

**原因**: SMILES 字符串格式不正确

**解决**:
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

### Q7: 找不到 Group 列

**错误信息**: `KeyError: 'Group'`

**解决**:
```r
# 在 R 中创建 Group 列
seurat_obj$Group <- seurat_obj$seurat_clusters  # 使用聚类结果
# 或
seurat_obj$Group <- seurat_obj$cell_type  # 使用细胞类型
```

---

## 附录：数据要求检查清单

### 转换前检查（R 端）

- [ ] 细胞数 ≥ 500
- [ ] 基因数已确定（将用于 `--gene_size` 参数）
- [ ] `Group` 列已创建
- [ ] 如使用药物结构：`SMILES` 和 `dose` 列已添加
- [ ] 已进行质量控制和标准化
- [ ] 数据使用正确的 slot (`data` 或 `counts`)

### 转换后验证（Python 端）

- [ ] `adata.n_obs`（细胞数）正确
- [ ] `adata.n_vars`（基因数）与 `--gene_size` 一致
- [ ] `adata.obs['Group']` 存在
- [ ] 如使用药物结构：`adata.obs['SMILES']` 和 `adata.obs['dose']` 存在
- [ ] SMILES 格式验证通过
- [ ] 无缺失值（NaN）
- [ ] 数据类型正确（ndarray 或 sparse）

### 训练前准备

- [ ] `--gene_size` 参数与 `adata.n_vars` 一致
- [ ] `--output_dim` 参数设置（通常等于 `gene_size`）
- [ ] 如使用药物结构：提供对照组数据路径
- [ ] `--batch_size` 根据细胞数合理设置

---

## 附录：快速参考

### R 端常用命令

```r
# 查看 Seurat 对象
print(seurat_obj)
print(paste("细胞数:", ncol(seurat_obj)))
print(paste("基因数:", nrow(seurat_obj)))
head(seurat_obj@meta.data)
colnames(seurat_obj@meta.data)

# 添加新列到元数据
seurat_obj$new_column <- values

# 筛选细胞
seurat_obj_subset <- subset(seurat_obj, condition == "treated")

# 筛选基因
seurat_obj <- seurat_obj[variable_genes, ]

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
print(f"细胞数: {adata.n_obs}")
print(f"基因数: {adata.n_vars}")
adata.obs.head()
adata.var.head()

# 筛选细胞
adata_subset = adata[adata.obs['condition'] == 'treated'].copy()

# 筛选基因
adata_subset = adata[:, adata.var['highly_variable']].copy()

# 保存
adata.write_h5ad("data.h5ad")
```

---

## 相关资源

- **Seurat 文档**: https://satijalab.org/seurat/
- **SeuratDisk 文档**: https://mojaveazure.github.io/seurat-disk/
- **AnnData 文档**: https://anndata.readthedocs.io/
- **Scanpy 文档**: https://scanpy.readthedocs.io/
- **SMILES 格式**: https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html

---

*文档最后更新：2025年*
