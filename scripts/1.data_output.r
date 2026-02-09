# --- 起手式 ---
print("Initiating High-End Visualization Pipeline...")
rm(list = ls())
setwd("E:/R/20231030_acu_colon/batch2_seurat/envs")
# 加载环境 (假设你有这个机制)
if(dir.exists("./")){
  lf <- list.files("./")
  for (file in lf) source(file)
}

# 设置工作目录
setwd("E:/Development/Squidiff")
dir.create("./data", showWarnings = FALSE)
setwd("./data")

library(Seurat)
library(SeuratDisk)

print("Loading data...")
# colon <- readRDS("E:/R/20231030_acu_colon/batch2_seurat/colon.rds")
# healthy_colon <- readRDS("E:/R/20231030_acu_colon/batch2_seurat/healthy_colon.rds")

fascia <- readRDS("E:/R/20231030_acu_colon/batch2_seurat/fascia.rds")
fascia$Group <- fascia$sample

table(fascia$cellclusters2)

#B cells              Dendritic cell           Endothelial cells
#27                         131                        1632
#Fibroblasts                Granulocytes Lymphatic endothelial cells
#5388                         362                         154
#Macrophages                  Mast cells              MyoFibroblasts
#5126                         402                        1283
#NK cells                 NKL T cells                Plasma cells
#189                        1270                          40
#T helper cells                   Telocytes
#232                        7260

table(fascia$Group)

#acupuncture   dysentery
#11185       12311

# ============================================================
# 步骤 1: 筛选 Mast cell 类型
# ============================================================
print("Step 1: Filtering Mast cells...")

mast_cells <- subset(fascia, subset = cellclusters2 == "Mast cells")
print(paste("Mast cells 数量:", ncol(mast_cells)))

# 查看分组情况
print("Group 分布:")
table(mast_cells$Group)

# ============================================================
# 步骤 2: 标准化数据
# ============================================================
print("Step 2: Normalizing data...")

mast_cells <- NormalizeData(
  mast_cells,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

mast_cells <- FindVariableFeatures(
  mast_cells,
  selection.method = "vst",
  nfeatures = 2000
)

# ============================================================
# 步骤 3: 计算差异基因 (top 500)
# ============================================================
print("Step 3: Finding marker genes (top 500)...")

# 设置 ident 为 Group
Idents(mast_cells) <- mast_cells$Group

# 计算差异基因
markers <- FindAllMarkers(
  mast_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# 按 avg_log2FC 排序，取 top 500
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 500, wt = avg_log2FC)

print(paste("Top 500 差异基因数量:", nrow(top_markers)))

# 获取唯一的基因列表
top_genes <- unique(top_markers$gene)
print(paste("唯一差异基因数量:", length(top_genes)))

# 保存差异基因列表
write.csv(top_markers, "mast_cell_top500_markers.csv", row.names = FALSE)
writeLines(top_genes, "mast_cell_top500_genes.txt")

# ============================================================
# 步骤 4: 筛选 top 500 差异基因
# ============================================================
print("Step 4: Subsetting to top 500 differential genes...")

# 如果 top_genes 超过 500，取前 500
if (length(top_genes) > 500) {
  top_genes <- top_genes[1:500]
}

# 筛选数据
mast_cells_top500 <- mast_cells[top_genes, ]
print(paste("筛选后基因数:", nrow(mast_cells_top500)))

# ============================================================
# 步骤 5: 进一步筛选为 200 个高变基因
# ============================================================
print("Step 5: Selecting top 200 highly variable genes...")

# 重新计算高变基因（在筛选后的 500 基因中）
mast_cells_top500 <- FindVariableFeatures(
  mast_cells_top500,
  selection.method = "vst",
  nfeatures = 500
)

# 获取高变基因
hv_genes <- VariableFeatures(mast_cells_top500)[1:200]
print(paste("选择的 200 个基因数量:", length(hv_genes)))

# 筛选为 200 个基因
mast_cells_200 <- mast_cells_top500[hv_genes, ]
print(paste("最终基因数:", nrow(mast_cells_200)))
print(paste("最终细胞数:", ncol(mast_cells_200)))

# 保存基因列表
writeLines(hv_genes, "mast_cell_200_genes.txt")

# ============================================================
# 步骤 6: 添加必需的 metadata
# ============================================================
print("Step 6: Adding required metadata...")

# 确保 Group 列存在（已存在）
# 可以根据需要添加其他 metadata

# ============================================================
# 步骤 7: 保存为 h5ad 格式
# ============================================================
print("Step 7: Saving to h5ad format...")

# 方法 1: 使用 SeuratDisk 保存为 h5ad
tryCatch({
  SaveH5Seurat(mast_cells_200, filename = "mast_cells_200genes.h5ad")
  print("✓ 数据已保存为: mast_cells_200genes.h5ad")
}, error = function(e) {
  # 方法 2: 如果 SaveH5Seurat 不可用，使用 Convert
  print("SaveH5Seurat 不可用，使用备用方法...")
  Save(mast_cells_200, filename = "mast_cells_200genes.h5seurat")
  Convert("mast_cells_200genes.h5seurat", dest = "h5ad", overwrite = TRUE)
  print("✓ 数据已保存为: mast_cells_200genes.h5ad")
})

Convert("mast_cells_200genes.h5ad.h5seurat", dest = "h5ad", overwrite = TRUE)
print("✓ 数据已保存为: mast_cells_200genes.h5ad")

# ============================================================
# 步骤 8: 数据验证
# ============================================================
print("Step 8: Data validation...")

print(paste("  - 细胞数:", ncol(mast_cells_200)))
print(paste("  - 基因数:", nrow(mast_cells_200)))
print(paste("  - Group 数量:", length(unique(mast_cells_200$Group))))
print("  - Group 分布:")
print(table(mast_cells_200$Group))

# 保存 metadata
metadata_df <- mast_cells_200@meta.data
write.csv(metadata_df, "mast_cells_200_metadata.csv", row.names = TRUE)

# ============================================================
# 完成
# ============================================================
print("=")
print("✓ 数据处理完成！")
print("  输出文件:")
print("    - mast_cells_200genes.h5ad (主数据文件)")
print("    - mast_cell_top500_markers.csv (差异基因列表)")
print("    - mast_cell_200_genes.txt (200个基因列表)")
print("    - mast_cells_200_metadata.csv (元数据)")
print("=")