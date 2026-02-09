import scanpy as sc
try:
    adata = sc.read_h5ad(r'E:\Development\Squidiff\data\mast_cells_200genes.h5ad')
    print(f'Data shape: {adata.shape}')
except Exception as e:
    print(e)
