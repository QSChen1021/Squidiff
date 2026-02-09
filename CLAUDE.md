# CLAUDE.md — Squidiff 开发指南

Squidiff 是一个基于扩散模型的单细胞转录组预测框架。这是一份 AI 辅助开发指南。

## 项目概述

**Squidiff** 用于预测细胞发育和对扰动的响应，支持：
- 药物处理后的单细胞转录组预测
- 细胞分化预测
- 基因扰动预测

## 技术栈

- **语言**: Python 3.8+
- **核心框架**: PyTorch (扩散模型实现)
- **数据处理**: AnnData (h5ad 格式)
- **构建工具**: setuptools
- **可选依赖**: PyQt5 (GUI), rich/click (CLI)

## 核心模块

| 模块 | 功能 |
|------|------|
| `diffusion.py` | 扩散模型核心实现 |
| `MLPModel.py` | MLP 神经网络架构 |
| `train_squidiff.py` | 训练脚本入口 |
| `sample_squidiff.py` | 采样/推理脚本 |
| `scrna_datasets.py` | 单细胞数据集处理 |
| `resample.py`, `respace.py` | 重采样与时间步调度 |
| `losses.py` | 损失函数定义 |

## 常用命令

### 训练模型
```bash
# 基础训练
python train_squidiff.py \
  --logger_path LOGGER_FIRE_NAME \
  --data_path YOUR_DATASET.h5ad \
  --resume_checkpoint ptNAME \
  --gene_size 500 \
  --output_dim 500

# 结合药物结构训练
python train_squidiff.py \
  --logger_path logger_files/logger_sciplex \
  --data_path datasets/sci_plex_train.h5ad \
  --resume_checkpoint sciplex_results \
  --use_drug_structure True \
  --gene_size 200 \
  --output_dim 200 \
  --control_data_path datasets/sci_plex_control.h5ad
```

### 模型推理
```python
import sample_squidiff

sampler = sample_squidiff.sampler(
    model_path='simu_results/model.pt',
    gene_size=100,
    output_dim=100,
    use_drug_structure=False
)

# 加载测试数据
import scanpy as sc
test_adata = sc.read_h5ad('datasets/test.h5ad')

# 获取编码并预测
z_sem = sampler.model.encoder(torch.tensor(test_adata.X).to('cuda'))
pred = sampler.pred(z_sem, gene_size=test_adata.shape[1])
```

## 开发规范

### 1. 代码组织
- **模型定义**: 放在 `Squidiff/` 模块内
- **训练脚本**: 根目录 `train_squidiff.py`
- **采样脚本**: 根目录 `sample_squidiff.py`
- **数据集**: `datasets/` 目录

### 2. 数据格式要求
输入 h5ad 文件需包含：
- 单细胞计数矩阵 (`.X`)
- 元数据 (`.obs`)
- (可选) 药物化合物信息

### 3. 调试原则
- 先检查数据维度和格式 (`gene_size`, `output_dim` 必须匹配)
- GPU 内存不足时减小 `batch_size` 或 `gene_size`
- 训练不稳定时检查学习率和扩散时间步设置

### 4. 扩展功能
添加新功能时遵循以下步骤：
1. **实现**: 在相应模块添加核心逻辑
2. **测试**: 使用小规模数据验证
3. **封装**: 整合到 `train_squidiff.py` 或 `sample_squidiff.py`
4. **文档**: 更新 README.md 参数说明

## 相关资源

- **复现仓库**: https://github.com/siyuh/Squidiff_reproducibility
- **论文引用**: He, S. et al. Squidiff: predicting cellular development and responses to perturbations using a diffusion model. Nat Methods (2025)

## 项目结构

```
Squidiff/
├── Squidiff/           # 核心包
│   ├── diffusion.py    # 扩散模型
│   ├── MLPModel.py     # 神经网络
│   └── ...
├── train_squidiff.py   # 训练入口
├── sample_squidiff.py  # 推理入口
├── datasets/           # 数据目录
└── README.md           # 项目说明
```
