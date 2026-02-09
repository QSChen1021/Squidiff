<div align="center">
  <img src="squidiff_logo.png" width="96" alt="Squidiff logo" />
  <h1>Squidiff</h1>
  <p><strong>Predicting cellular development and perturbation responses with diffusion models</strong></p>
  <p><strong>基于扩散模型预测细胞发育与扰动响应</strong></p>
</div>

[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/)
[![PyTorch 2.0+](https://img.shields.io/badge/PyTorch-2.0%2B-ee4c2c.svg)](https://pytorch.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## 1. Overview / 项目概览

**EN:** Squidiff is a diffusion-model framework for single-cell transcriptomics.  
**中文：** Squidiff 是一个面向单细胞转录组任务的扩散模型框架。  

**EN:** It learns transitions from baseline cellular states to perturbed states and generates predicted expression profiles.  
**中文：** 它学习细胞从基线状态到扰动状态的转换过程，并生成预测表达谱。  

## 2. What You Can Do / 能做什么

- **EN:** Predict transcriptomic responses to drugs and gene perturbations.  
  **中文：** 预测药物和基因扰动后的转录组响应。  
- **EN:** Model cell-state trajectories in development and differentiation tasks.  
  **中文：** 建模细胞发育与分化过程中的状态轨迹。  
- **EN:** Integrate optional molecular structure features (SMILES + dose).  
  **中文：** 支持可选的分子结构特征（SMILES + dose）融合。  

## 3. Installation / 安装

### 3.1 PyPI

```bash
pip install Squidiff
```

### 3.2 Local Development (Recommended) / 本地开发（推荐）

```bash
git clone https://github.com/siyuh/Squidiff.git
cd Squidiff
```

```bash
# Windows: create venv with Python 3.11
py -3.11 -m venv .venv

# Linux/macOS: create venv with Python 3.11
python3.11 -m venv .venv
```

```bash
# Windows: activate
.venv\Scripts\activate

# Linux/macOS: activate
source .venv/bin/activate
```

```bash
pip install -e ".[dev]"
```

**EN:** If you prefer dependency lock by file, install via `requirements.txt`.  
**中文：** 如果你偏好按文件固定依赖，也可使用 `requirements.txt` 安装。  

```bash
pip install -r requirements.txt
```

## 4. Data Requirements / 数据要求

**EN:** Input format is `h5ad` (`AnnData`).  
**中文：** 输入格式为 `h5ad`（`AnnData`）。  

| Field | Required | EN | 中文 |
|---|---|---|---|
| `adata.X` | Yes | Cell-by-gene expression matrix | 细胞 x 基因表达矩阵 |
| `adata.obs["Group"]` | Yes | Group/condition label | 分组或条件标签 |
| `adata.obs["SMILES"]` | Conditional | Required when `use_drug_structure=True` | 当 `use_drug_structure=True` 时必需 |
| `adata.obs["dose"]` | Conditional | Required when `use_drug_structure=True` | 当 `use_drug_structure=True` 时必需 |

## 5. Quick Start / 快速开始

### 5.1 Check Input Shape / 检查输入维度

```bash
python scripts/check_shape.py --data_path path/to/train.h5ad
```

**EN:** Use reported gene count for both `--gene_size` and `--output_dim`.  
**中文：** 将输出的基因数同时用于 `--gene_size` 与 `--output_dim`。  

### 5.2 Train Baseline Model / 训练基础模型

```bash
python train_squidiff.py \
  --logger_path "./results/baseline" \
  --data_path "path/to/train.h5ad" \
  --resume_checkpoint "./checkpoints/baseline" \
  --gene_size 159 \
  --output_dim 159
```

### 5.3 Train with Drug Structure / 使用药物结构训练

```bash
python train_squidiff.py \
  --logger_path "./results/drug" \
  --data_path "path/to/train.h5ad" \
  --control_data_path "path/to/control.h5ad" \
  --resume_checkpoint "./checkpoints/drug" \
  --use_drug_structure True \
  --gene_size 159 \
  --output_dim 159
```

### 5.4 Inference / 推理

```python
import torch
import scanpy as sc
import sample_squidiff

sampler = sample_squidiff.sampler(
    model_path="checkpoints/baseline/model.pt",
    gene_size=159,
    output_dim=159,
    use_drug_structure=False,
)

test_adata = sc.read_h5ad("path/to/test.h5ad")
z_sem = sampler.model.encoder(torch.tensor(test_adata.X).to("cuda"))
pred = sampler.pred(z_sem, gene_size=test_adata.shape[1])
```

## 6. Key Training Arguments / 核心训练参数

| Argument | EN | 中文 |
|---|---|---|
| `--logger_path` | Output directory for logs/checkpoints | 日志与检查点输出目录 |
| `--data_path` | Training dataset path (`h5ad`) | 训练数据路径（`h5ad`） |
| `--resume_checkpoint` | Resume/save checkpoint location | 恢复/保存检查点位置 |
| `--gene_size` | Input gene dimension | 输入基因维度 |
| `--output_dim` | Output dimension (usually same as `gene_size`) | 输出维度（通常与 `gene_size` 相同） |
| `--use_drug_structure` | Enable drug-structure conditioning | 开启药物结构条件建模 |
| `--control_data_path` | Control dataset for drug-structure mode | 药物结构模式下的对照组数据路径 |
| `--batch_size` | Training batch size | 训练批大小 |
| `--lr` | Learning rate | 学习率 |

## 7. Documentation Map / 文档导航

- **EN:** Deployment + environment + training operations: `docs/部署文档.md`  
  **中文：** 部署、环境、训练运行文档：`docs/部署文档.md`  
- **EN:** Seurat to h5ad conversion and data validation: `docs/seurat转换指南.md`  
  **中文：** Seurat 转换与数据校验文档：`docs/seurat转换指南.md`  
- **EN:** Troubleshooting checklist and common errors: `docs/避坑指南.md`  
  **中文：** 排错清单与常见错误：`docs/避坑指南.md`  

## 8. Quality Check / 质量检查

```bash
python -m ruff check .
```

## 9. Reproducibility / 复现资源

**EN:** Additional end-to-end resources are available at:  
**中文：** 更多端到端复现资源请见：  

- https://github.com/siyuh/Squidiff_reproducibility

## 10. Citation / 引用

```bibtex
@article{he2025squidiff,
  title={Squidiff: predicting cellular development and responses to perturbations using a diffusion model},
  author={He, Siyu and Zhu, Yitan and Tavakol, Diana N and others},
  journal={Nature Methods},
  year={2025},
  doi={10.1038/s41592-025-02877-y}
}
```

## 11. Contact / 联系方式

- Siyu He: siyuhe@stanford.edu
- GitHub Issues: https://github.com/siyuh/Squidiff/issues

## 12. License / 许可证

**EN:** Released under the MIT License.  
**中文：** 本项目采用 MIT 许可证。  
