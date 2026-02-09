# 脚本验证报告 - 2025-02-09

> 系统性检查 `E:\Development\Squidiff\scripts\` 目录下的所有 Python 脚本

---

## 检查概述

| 检查项 | 状态 |
|--------|------|
| Ruff 代码质量检查 | ✅ 全部通过 |
| Python 语法检查 | ✅ 全部通过 |
| 导入依赖检查 | ⚠️ 部分缺失 |

---

## 脚本列表

### 1. `scripts/3.sample_mast_cells.py`

**功能**: Squidiff 采样/推理脚本 - Mast Cell 数据

**状态**: ✅ 代码质量合格

**修复问题**:
- 移除未使用的 `numpy` 导入 (F401)
- 修复无占位符的 f-string (F541)

**依赖要求**:
```
torch>=2.0.0      [MISSING] - pip install torch
scanpy>=1.9.0     [MISSING] - pip install scanpy
```

**可选依赖**:
```
scikit-learn      (用于 R² 计算)
```

---

### 2. `scripts/4.validate_data.py`

**功能**: 数据验证脚本 - 验证 h5ad 数据是否符合 Squidiff 要求

**状态**: ✅ 代码质量合格

**修复问题**:
- 移除重复的 docstring
- 修复 4 个无占位符的 f-string (F541)
- 修复 bare except 为 `except Exception:` (E722)

**依赖要求**:
```
scanpy>=1.9.0     [MISSING] - pip install scanpy
```

**可选依赖**:
```
rdkit>=2023.3.1   [INFO] 未安装 (可选，用于 SMILES 验证)
```

---

### 3. `scripts/5.convert_h5seurat_to_h5ad.py`

**功能**: 将 Seurat 导出的 h5seurat 文件转换为 h5ad 格式

**状态**: ✅ 代码质量合格

**修复问题**:
- 移除未使用的 `numpy` 导入 (F401)

**依赖要求**:
```
h5py>=3.8.0       ✅ [INSTALLED]
pandas            ✅ [INSTALLED]
anndata>=0.8.0    [MISSING] - pip install anndata
```

**可选依赖**:
```
seurat-disk       [OPTIONAL] - 备用转换方法
```

---

## 依赖安装命令

```bash
# 安装所有必需依赖
pip install torch>=2.0.0 scanpy>=1.9.0 anndata>=0.8.0 h5py>=3.8.0 pandas scikit-learn

# 安装可选依赖 (药物结构处理)
pip install rdkit>=2023.3.1
# 或使用 uv
uv pip install rdkit>=2023.3.1
```

---

## 脚本使用流程

### 完整流水线

```
1. data_output.r (R) → 生成 h5seurat 文件
   ↓
2. convert_h5seurat_to_h5ad.py → 转换为 h5ad
   ↓
3. validate_data.py → 验证数据格式
   ↓
4. train_mast_cells.sh/.bat → 训练模型
   ↓
5. sample_mast_cells.py → 推理预测
```

### 快速命令参考

```bash
# 1. 格式转换 (Python)
python scripts/5.convert_h5seurat_to_h5ad.py data/mast_cells_200genes.h5seurat

# 2. 数据验证
python scripts/4.validate_data.py --data_path data/mast_cells_200genes.h5ad

# 3. 模型训练 (Windows)
scripts\2.train_mast_cells.bat

# 4. 模型推理
python scripts/3.sample_mast_cells.py \
  --model_path logs/mast_cells_<timestamp>/model.pt \
  --data_path data/mast_cells_200genes.h5ad \
  --output_dir outputs
```

---

## Ruff 修复详情

### 修复前统计
- 错误数: 8 个
- 可自动修复: 7 个

### 修复后统计
- 错误数: 0 个 ✅

### 错误类型分布
| 错误代码 | 数量 | 说明 |
|----------|------|------|
| F401 | 2 | 未使用的导入 |
| F541 | 5 | 无占位符的 f-string |
| E722 | 1 | 裸 except |

---

## 下一步建议

1. **安装缺失依赖**: 运行上述 pip 命令安装所有必需包
2. **测试转换脚本**: 用实际的 h5seurat 文件测试转换功能
3. **端到端测试**: 运行完整的 R → Python → 训练 → 推理流程

---

*检查时间: 2025-02-09*
*检查工具: ruff, py_compile, import check*
