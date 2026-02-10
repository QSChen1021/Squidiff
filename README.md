<div align="center">
  <img src="squidiff_logo.png" width="80" />
  <h1>Squidiff</h1>
  <p>
    <strong>Predicting cellular development and responses to perturbations using a diffusion model</strong>
  </p>
  <p>
    <strong>基于扩散模型的单细胞转录组预测框架 —— 预测细胞发育和对扰动的响应</strong>
  </p>
</div>

---

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![PyTorch](https://img.shields.io/badge/PyTorch-1.10+-ee4c2c.svg)](https://pytorch.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Overview / 概述

**English** | [中文](#中文)

---

### What is Squidiff? / Squidiff 是什么？

Squidiff is a **diffusion model-based generative framework** designed to predict single-cell transcriptomic changes across diverse cell types in response to a wide range of perturbations.

**Squidiff** 是一个**基于扩散模型的生成式框架**，用于预测多种细胞类型在不同扰动条件下的单细胞转录组变化。

---

### Key Features / 核心功能

| Feature | 功能 | Description |
|---------|------|-------------|
| 🧪 **Drug Response Prediction** | 药物响应预测 | Predict transcriptomic changes after drug treatments |
| 🧬 **Cell Differentiation** | 细胞分化预测 | Model cell development trajectories and fate decisions |
| 🔬 **Gene Perturbation** | 基因扰动预测 | Simulate effects of gene knockouts/overexpression |
| 💊 **Drug Structure Integration** | 药物结构整合 | Incorporate molecular structures (SMILES) for better predictions |

---

### How It Works / 工作原理

```
┌─────────────────────────────────────────────────────────────────┐
│                    Training Phase / 训练阶段                     │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│   Input Data / 输入数据                                          │
│   ┌──────────────┐         ┌──────────────┐                    │
│   │   Before     │         │    After      │                    │
│   │  Perturbation│  ──────> │  Perturbation│                    │
│   │  (e.g. untreated)     │  (e.g. drug-treated)            │
│   └──────────────┘         └──────────────┘                    │
│                                                                  │
│   ↓                                                              │
│   Diffusion Model Learns / 扩散模型学习                          │
│   "How to model the transition between cellular states"         │
│   "如何建模细胞状态之间的转换"                                   │
│                                                                  │
│   ↓                                                              │
│   Trained Model / 训练好的模型                                   │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                  Inference Phase / 推理阶段                      │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│   New Cell + Perturbation / 新细胞 + 扰动                       │
│   ┌──────────────┐         ┌──────────────┐                    │
│   │   Cell State │   +     │  Drug / Gene  │                    │
│   │              │         │  Perturbation │                    │
│   └──────────────┘         └──────────────┘                    │
│                                                                  │
│   ↓                                                              │
│   Squidiff Model / Squidiff 模型                                │
│                                                                  │
│   ↓                                                              │
│   Predicted Transcriptomic Response / 预测的转录组响应           │
└─────────────────────────────────────────────────────────────────┘
```

---

### Model Architecture / 模型架构

```
Squidiff Architecture:

Input (Gene Expression)       ──┐
                                    │
Timestep Embedding             ──┼──> MLP Encoder
                                    │        │
Optional: Drug Structure       ──┘        │
(Optional: Cell Group)                      │
                                             ▼
                                      MLP Blocks (×N)
                                      with LayerNorm
                                             │
                                             ▼
                                      Output Layer
                                             │
                                             ▼
                              Predicted Gene Expression
```

**Key Components / 核心组件:**

- **MLP Encoder**: Encodes cell state into latent representation
- **Timestep Embedding**: Diffusion timestep conditioning
- **Drug Encoder** (optional): Encodes SMILES + dose information
- **MLP Decoder**: Generates predicted gene expression

---

## Installation / 安装

```bash
pip install Squidiff
```

### Dependencies / 依赖项

```
python >= 3.8
torch >= 1.10
scanpy
anndata
h5py
numpy
pandas
```

For drug structure integration / 药物结构整合:
```
rdkit
```

---

## Quick Start / 快速开始

### 1. Prepare Your Data / 准备数据

Prepare an `h5ad` file with / 准备包含以下内容的 `h5ad` 文件:

- **Single-cell count matrix** (`adata.X`): Cells × Genes
- **Metadata** (`adata.obs`): Cell annotations (must include `Group` column)
- **Optional** (`adata.obs`): `SMILES` (drug structure), `dose`

### 2. Training / 训练

#### Basic Training (No drug structure) / 基础训练（无药物结构）

```bash
python train_squidiff.py \
  --logger_path "./results" \
  --data_path "data/mast_cells.h5ad" \
  --resume_checkpoint "./checkpoints" \
  --gene_size 159 \
  --output_dim 159
```

**Important Parameters / 重要参数:**

| Parameter | Description | 说明 |
|-----------|-------------|------|
| `--gene_size` | Number of genes in dataset | 数据集中的基因数量 |
| `--output_dim` | Output dimension (should = gene_size) | 输出维度（应等于基因数） |
| `--logger_path` | Directory for logs | 日志保存目录 |
| `--resume_checkpoint` | Directory for model checkpoints | 模型检查点保存目录 |

#### Training with Drug Structure / 结合药物结构训练

```bash
python train_squidiff.py \
  --logger_path "./logs" \
  --data_path "datasets/sciplex_train.h5ad" \
  --resume_checkpoint "./checkpoints" \
  --use_drug_structure True \
  --gene_size 200 \
  --output_dim 200 \
  --control_data_path "datasets/sciplex_control.h5ad"
```

### 3. Sampling / Inference / 采样/推理

```python
import sample_squidiff
import scanpy as sc
import torch

# Initialize sampler / 初始化采样器
sampler = sample_squidiff.sampler(
    model_path='checkpoints/model.pt',
    gene_size=159,
    output_dim=159,
    use_drug_structure=False
)

# Load test data / 加载测试数据
test_adata = sc.read_h5ad('datasets/test.h5ad')

# Get latent encoding / 获取编码
z_sem = sampler.model.encoder(
    torch.tensor(test_adata.X).to('cuda')
)

# Predict / 预测
predicted_expression = sampler.pred(
    z_sem,
    gene_size=test_adata.shape[1]
)
```

---

## Project Structure / 项目结构

```
Squidiff/
├── Squidiff/                    # Core package / 核心包
│   ├── diffusion.py             # Diffusion model / 扩散模型
│   ├── MLPModel.py              # Neural network architecture / 神经网络架构
│   ├── scrna_datasets.py        # Data loading / 数据加载
│   ├── train_util.py            # Training utilities / 训练工具
│   ├── resample.py              # Resampling functions / 重采样
│   ├── respace.py               # Timestep spacing / 时间步调度
│   └── losses.py                # Loss functions / 损失函数
│
├── train_squidiff.py            # Training script / 训练脚本
├── sample_squidiff.py           # Sampling script / 采样脚本
│
├── scripts/                     # Utility scripts / 工具脚本
│   ├── 0.pipeline_summary.md    # Pipeline documentation / 流程文档
│   ├── 1.data_output.r          # R data export / R 数据导出
│   ├── 2.train_mast_cells.sh    # Training example / 训练示例
│   ├── 3.sample_mast_cells.py   # Sampling example / 采样示例
│   ├── 4.validate_data.py       # Data validation / 数据验证
│   └── 5.convert_h5seurat_to_h5ad.py  # Format conversion / 格式转换
│
├── data/                        # Data directory / 数据目录
├── logs/                        # Training logs / 训练日志
├── outputs/                     # Inference results / 推理结果
└── README.md
```

---

## Use Cases / 应用场景

### Scenario 1: Drug Discovery / 场景 1：药物发现

```
Question: What will this new drug do to cells?
问题: 这个新药对细胞会有什么影响？

Traditional: Wet lab experiments → weeks/months
传统方法: 湿实验 → 数周到数月

Squidiff: Input drug structure → seconds
Squidiff: 输入药物结构 → 秒级预测
```

### Scenario 2: Cell Fate Prediction / 场景 2：细胞命运预测

```
Question: How will stem cells differentiate?
问题: 干细胞会如何分化？

Squidiff: Predict intermediate states along developmental trajectories
Squidiff: 预测发育轨迹上的中间状态
```

### Scenario 3: Gene Perturbation / 场景 3：基因扰动

```
Question: What happens if I knock out Gene X?
问题: 如果敲除基因 X 会发生什么？

Squidiff: Predict transcriptomic consequences of genetic perturbations
Squidiff: 预测基因扰动的转录组后果
```

---

## Reproducibility / 复现性

For complete data preparation, model usage, and downstream analysis examples, please visit:

完整的数据准备、模型使用和下游分析示例，请访问：

**https://github.com/siyuh/Squidiff_reproducibility**

---

## Citation / 引用

If you use Squidiff in your research, please cite:

如果您在研究中使用了 Squidiff，请引用：

```bibtex
@article{he2025squidiff,
  title={Squidiff: predicting cellular development and responses to perturbations using a diffusion model},
  author={He, Siyu and Zhu, Yitan and Tavakol, Diana N and others},
  journal={Nature Methods},
  year={2025},
  doi={10.1038/s41592-025-02877-y}
}

@article{nature2025squidiff,
  title={Predicting cellular responses with conditional diffusion models},
  journal={Nature Methods},
  year={2025},
  doi={10.1038/s41592-025-02878-x}
}
```

---

## Contact / 联系方式

**Questions? / 有问题?**

- **Siyu He** - siyuhe@stanford.edu
- **GitHub Issues** - [Create an issue](https://github.com/siyuh/Squidiff/issues)

---

## License / 许可证

This project is licensed under the MIT License - see the LICENSE file for details.

本项目采用 MIT 许可证 - 详见 LICENSE 文件。

---

<div align="center">
  <p>Built with ❤️ for single-cell research</p>
  <p>为单细胞研究而构建</p>
</div>
