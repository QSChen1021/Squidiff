#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
============================================================
数据验证脚本 - 验证 h5ad 数据是否符合 Squidiff 要求
============================================================
"""

import os
import sys
import argparse
import scanpy as sc

# 可选导入 rdkit
try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("警告: rdkit 未安装，将跳过 SMILES 验证")


def validate_data(data_path, use_drug_structure=False):
    """
    验证 AnnData 对象是否符合 Squidiff 要求

    参数:
        data_path: h5ad 文件路径
        use_drug_structure: 是否使用药物结构

    返回:
        is_valid: bool
        issues: list of str
    """
    issues = []

    # 检查是否使用药物结构但 rdkit 不可用
    if use_drug_structure and not RDKIT_AVAILABLE:
        issues.append("❌ 使用药物结构需要安装 rdkit: conda install -c conda-forge rdkit")
        return False, issues
    """
    验证 AnnData 对象是否符合 Squidiff 要求

    参数:
        data_path: h5ad 文件路径
        use_drug_structure: 是否使用药物结构

    返回:
        is_valid: bool
        issues: list of str
    """
    issues = []

    print("=" * 60)
    print("数据验证脚本")
    print("=" * 60)
    print(f"文件路径: {data_path}")
    print()

    # 检查文件是否存在
    if not os.path.exists(data_path):
        print(f"❌ 文件不存在: {data_path}")
        return False, ["文件不存在"]

    # 读取数据
    try:
        adata = sc.read_h5ad(data_path)
    except Exception as e:
        print(f"❌ 读取文件失败: {e}")
        return False, [f"读取文件失败: {e}"]

    # 检查基本结构
    print("基本结构:")
    print(f"  ✓ 细胞数: {adata.n_obs}")
    print(f"  ✓ 基因数: {adata.n_vars}")
    print(f"  ✓ 矩阵类型: {type(adata.X)}")

    # 检查维度
    n_cells, n_genes = adata.n_obs, adata.n_vars

    if n_cells < 500:
        issues.append(f"⚠️  细胞数较少 ({n_cells})，建议至少 500 个")

    if n_genes < 50:
        issues.append(f"⚠️  基因数过少 ({n_genes})，建议至少 100 个")

    # 检查必需列
    print("\n必需列检查:")
    if 'Group' not in adata.obs.columns:
        issues.append("❌ 缺少必需列: Group")
    else:
        print(f"  ✓ Group 列存在，{adata.obs['Group'].nunique()} 个唯一值")
        print(f"    分布:")
        for group, count in adata.obs['Group'].value_counts().items():
            print(f"      - {group}: {count}")

    # 检查药物结构相关列
    if use_drug_structure:
        print("\n药物结构检查:")
        if 'SMILES' not in adata.obs.columns:
            issues.append("❌ 使用药物结构时缺少必需列: SMILES")
        else:
            # 验证 SMILES 格式
            if RDKIT_AVAILABLE:
                valid_count = 0
                for smiles in adata.obs['SMILES']:
                    try:
                        if Chem.MolFromSmiles(smiles) is not None:
                            valid_count += 1
                    except:
                        pass

                if valid_count < len(adata):
                    issues.append(f"⚠️  有 {len(adata) - valid_count} 个无效的 SMILES")
                else:
                    print(f"  ✓ SMILES 格式正确 ({valid_count}/{len(adata)})")
            else:
                print("  ⚠️  rdkit 未安装，跳过 SMILES 格式验证")

        if 'dose' not in adata.obs.columns:
            issues.append("❌ 使用药物结构时缺少必需列: dose")
        else:
            dose_range = (adata.obs['dose'].min(), adata.obs['dose'].max())
            print(f"  ✓ dose 列存在，范围: {dose_range[0]:.2f} - {dose_range[1]:.2f}")

    # 检查数据类型
    print("\n数据类型检查:")
    import numpy as np
    from scipy.sparse import issparse

    if isinstance(adata.X, np.ndarray):
        print(f"  ✓ 密集矩阵格式")
    elif issparse(adata.X):
        print(f"  ✓ 稀疏矩阵格式")
    else:
        issues.append(f"⚠️  未知的数据类型: {type(adata.X)}")

    # 检查缺失值
    if hasattr(adata.X, 'data'):
        data_values = adata.X.data
    else:
        data_values = adata.X

    if np.isnan(data_values).any():
        issues.append("⚠️  存在缺失值 (NaN)")
    else:
        print(f"  ✓ 无缺失值")

    # 检查数据值范围
    print("\n数据值范围:")
    if hasattr(adata.X, 'data'):
        min_val = adata.X.data.min()
        max_val = adata.X.data.max()
    else:
        min_val = adata.X.min()
        max_val = adata.X.max()

    print(f"  - 最小值: {min_val:.4f}")
    print(f"  - 最大值: {max_val:.4f}")

    # 总结
    print("\n" + "=" * 60)
    if len(issues) == 0:
        print("✓ 数据验证通过！")
        print("=" * 60)
        return True, []
    else:
        print("发现问题:")
        for issue in issues:
            print(issue)
        print("=" * 60)
        return False, issues


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='验证 Squidiff 数据格式')
    parser.add_argument('--data_path', type=str, required=True,
                       help='h5ad 文件路径')
    parser.add_argument('--use_drug_structure', action='store_true',
                       help='是否使用药物结构')

    args = parser.parse_args()

    is_valid, issues = validate_data(
        args.data_path,
        use_drug_structure=args.use_drug_structure
    )

    sys.exit(0 if is_valid else 1)


if __name__ == "__main__":
    main()
