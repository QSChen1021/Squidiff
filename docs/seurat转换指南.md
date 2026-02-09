# Seurat 转换指南（R -> h5ad）

> 目标：将 Seurat 对象转换为 Squidiff 可直接使用的 `h5ad`，并完成最小验证。

## 1. 必要数据要求

转换后数据需满足：

- `adata.X` 存在，且维度为 `(n_cells, n_genes)`
- `adata.obs["Group"]` 存在
- 若使用药物结构：`adata.obs["SMILES"]` 与 `adata.obs["dose"]` 必须存在

训练参数约束：

- `--gene_size == adata.n_vars`
- `--output_dim == adata.n_vars`

## 2. 推荐转换路径（SeuratDisk）

### 2.1 R 端导出

```r
library(Seurat)
library(SeuratDisk)

# 确保 Group 存在
if (!"Group" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$Group <- seurat_obj$seurat_clusters
}

# 推荐流程：先 h5seurat，再转 h5ad
SaveH5Seurat(seurat_obj, filename = "sample.h5seurat", overwrite = TRUE)
Convert("sample.h5seurat", dest = "h5ad", overwrite = TRUE)
```

### 2.2 Python 端验证

```python
import scanpy as sc

adata = sc.read_h5ad("sample.h5ad")
print("cells:", adata.n_obs)
print("genes:", adata.n_vars)
print("obs columns:", adata.obs.columns.tolist())

assert "Group" in adata.obs.columns, "missing Group"
```

## 3. 药物结构场景补充

当训练使用 `--use_drug_structure True` 时：

1. 处理组数据中必须有 `SMILES` 和 `dose`。  
2. 需要提供 `--control_data_path`。  
3. 对照组与处理组的基因维度必须一致。  

SMILES 简单校验示例：

```python
from rdkit import Chem

def valid_smiles(s):
    return Chem.MolFromSmiles(s) is not None
```

## 4. 常见问题

### 4.1 `KeyError: 'Group'`

- 在 R 中补齐 `Group` 列后重新导出。

### 4.2 维度不匹配

- 用 `adata.n_vars` 回填训练命令里的 `--gene_size` 与 `--output_dim`。

### 4.3 `SMILES` 无效

- 用 RDKit 先批量过滤非法 SMILES，再训练。

## 5. 转换后检查清单

- [ ] `adata.n_obs > 0`
- [ ] `adata.n_vars > 0`
- [ ] `Group` 列存在
- [ ] 若使用药物结构：`SMILES` 与 `dose` 列存在
- [ ] `--gene_size` 与 `--output_dim` 已和 `adata.n_vars` 对齐
- [ ] 训练前可通过 `python scripts/check_shape.py --data_path ...` 复核
