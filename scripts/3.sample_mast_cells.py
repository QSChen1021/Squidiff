#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
============================================================
Squidiff 采样/推理脚本 - Mast Cell 数据
============================================================

用途: 使用训练好的 Mast Cell 模型进行预测
数据: data/mast_cells_200genes.h5ad
"""

import os
import sys
import argparse
import numpy as np
import torch
import scanpy as sc

# 添加项目路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sample_squidiff


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='Squidiff 采样脚本')

    # 模型参数
    parser.add_argument('--model_path', type=str, required=True,
                       help='模型文件路径')
    parser.add_argument('--gene_size', type=int, default=200,
                       help='基因数 (默认: 200)')
    parser.add_argument('--output_dim', type=int, default=200,
                       help='输出维度 (默认: 200)')
    parser.add_argument('--use_drug_structure', type=bool, default=False,
                       help='是否使用药物结构 (默认: False)')

    # 数据参数
    parser.add_argument('--data_path', type=str, required=True,
                       help='测试数据路径 (h5ad)')
    parser.add_argument('--output_dir', type=str, default='outputs',
                       help='输出目录 (默认: outputs)')

    # 采样参数
    parser.add_argument('--num_samples', type=int, default=100,
                       help='生成样本数量 (默认: 100)')
    parser.add_argument('--batch_size', type=int, default=16,
                       help='采样批次大小 (默认: 16)')

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_args()

    print("=" * 60)
    print("Squidiff 采样脚本 - Mast Cell")
    print("=" * 60)
    print(f"模型路径:     {args.model_path}")
    print(f"数据路径:     {args.data_path}")
    print(f"基因数:       {args.gene_size}")
    print(f"输出维度:     {args.output_dim}")
    print(f"生成样本数:   {args.num_samples}")
    print("=" * 60)
    print()

    # 检查模型文件
    if not os.path.exists(args.model_path):
        print(f"错误: 模型文件不存在: {args.model_path}")
        sys.exit(1)

    # 检查数据文件
    if not os.path.exists(args.data_path):
        print(f"错误: 数据文件不存在: {args.data_path}")
        sys.exit(1)

    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)

    # ============================================================
    # 步骤 1: 加载模型
    # ============================================================
    print("步骤 1: 加载模型...")

    sampler = sample_squidiff.sampler(
        model_path=args.model_path,
        gene_size=args.gene_size,
        output_dim=args.output_dim,
        use_drug_structure=args.use_drug_structure
    )
    print("✓ 模型加载完成")
    print()

    # ============================================================
    # 步骤 2: 加载测试数据
    # ============================================================
    print("步骤 2: 加载测试数据...")

    test_adata = sc.read_h5ad(args.data_path)
    print(f"  - 细胞数: {test_adata.n_obs}")
    print(f"  - 基因数: {test_adata.n_vars}")
    print(f"  - Group 数: {test_adata.obs['Group'].nunique()}")
    print(f"  - Group 分布:")
    print(test_adata.obs['Group'].value_counts())
    print("✓ 数据加载完成")
    print()

    # ============================================================
    # 步骤 3: 获取编码表示
    # ============================================================
    print("步骤 3: 获取编码表示...")

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"  - 设备: {device}")

    # 转换数据为 tensor
    if hasattr(test_adata.X, 'toarray'):
        x_data = torch.tensor(test_adata.X.toarray(), dtype=torch.float32).to(device)
    else:
        x_data = torch.tensor(test_adata.X, dtype=torch.float32).to(device)

    # 获取编码
    z_sem = sampler.model.encoder(x_data)
    print(f"  - 编码形状: {z_sem.shape}")
    print("✓ 编码完成")
    print()

    # ============================================================
    # 步骤 4: 生成预测
    # ============================================================
    print("步骤 4: 生成预测...")

    pred = sampler.pred(z_sem, gene_size=args.gene_size)
    print(f"  - 预测形状: {pred.shape}")
    print("✓ 预测完成")
    print()

    # ============================================================
    # 步骤 5: 保存结果
    # ============================================================
    print("步骤 5: 保存结果...")

    # 将预测结果添加到 AnnData
    test_adata.obsm['predicted'] = pred.detach().cpu().numpy()

    # 保存为 h5ad
    output_path = os.path.join(args.output_dir, 'mast_cells_predicted.h5ad')
    test_adata.write_h5ad(output_path)
    print(f"✓ 预测结果已保存: {output_path}")
    print()

    # ============================================================
    # 步骤 6: 计算评估指标
    # ============================================================
    print("步骤 6: 计算评估指标...")

    # 计算 MSE
    mse = torch.mean((x_data - pred[:x_data.shape[0]]) ** 2).item()
    print(f"  - MSE: {mse:.6f}")

    # 计算 R²
    from sklearn.metrics import r2_score
    x_np = x_data.cpu().numpy()
    pred_np = pred[:x_data.shape[0]].detach().cpu().numpy()
    r2 = r2_score(x_np.flatten(), pred_np.flatten())
    print(f"  - R²: {r2:.6f}")
    print("✓ 评估完成")
    print()

    # ============================================================
    # 完成
    # ============================================================
    print("=" * 60)
    print("✓ 采样完成！")
    print(f"输出文件: {output_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()
