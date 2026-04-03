# ===== 细胞注释分析 =====
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# 1. 加载并行包
library(future)

# 2. 设置并行策略
plan("multisession", workers = 4)

# 3. 增大内存限制
options(future.globals.maxSize = 50 * 1024^3)

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

setwd("/root/autodl-tmp/Single _cell")
load("seurat_object.rdata")

cat("细胞数:", ncol(seurat_object), "\n")
cat("Cluster数:", length(levels(seurat_object)), "\n")
cat("Clusters:", levels(seurat_object@active.ident), "\n")

# 创建graph目录
if (!dir.exists("graph")) dir.create("graph")

# ===== 方法一：手动注释 - 基于marker基因 =====
cat("\n===== 1. 绘制marker基因表达 =====\n")

# 定义的marker基因（来自原始脚本）
markers <- c("ALPL","RUNX2","IBSP",  # Osteoblastic OS cells
             "LYZ","CD68",         # Myeloid cells
             "ACP5","CTSK",        # OCs (Osteoclasts)
             "COL1A1","FAP","VIM", # CAFs
             "CD2","CD3D","CD3E","CD3G","GNLY","NKG7","KLRD1","KLRB1", # NK/T cells
             "EGFL7","PLVAP",     # Endothelial cells
             "MS4A1","CD79A",     # B cells
             "IGHG1","MZB1")      # Plasma cells

# 筛选存在的marker
available_markers <- markers[markers %in% rownames(seurat_object)]
cat("可用marker:", length(available_markers), "/", length(markers), "\n")

# FeaturePlot
p1 <- FeaturePlot(seurat_object, features = available_markers, ncol = 6)
ggsave("graph/Markers1.tiff", width = 16, height = 12, dpi = 300, plot = p1)

# DotPlot
p2 <- DotPlot(seurat_object, features = available_markers) + theme_bw() + RotatedAxis()
ggsave("graph/Markers2.tiff", width = 8, height = 6, dpi = 300, plot = p2)

# VlnPlot
p3 <- VlnPlot(seurat_object, features = available_markers, stack = TRUE, flip = TRUE) + NoLegend()
ggsave("graph/Markers3.tiff", width = 8, height = 12, dpi = 300, plot = p3)

# 查看cluster分布
p <- DimPlot(seurat_object, label = TRUE, reduction = "umap")
ggsave("graph/Cluster_Labels.tiff", width = 8, height = 6, dpi = 300, plot = p)

# ===== 方法二：自动找marker =====
cat("\n===== 2. 寻找各cluster的marker基因 =====\n")

allmarkers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(allmarkers, "allmarkers.csv")
cat("找到", nrow(allmarkers), "个marker基因\n")

# 每个cluster的前10个marker
Top10 <- allmarkers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(Top10, "Top10_coarse.csv")

# 显示每个cluster的top3 marker
cat("\n各cluster的top3 marker:\n")
for (c in unique(allmarkers$cluster)) {
  top3 <- allmarkers %>% filter(cluster == c) %>% slice_max(n = 3, order_by = avg_log2FC)
  cat("Cluster", c, ":", paste(top3$gene, collapse = ", "), "\n")
}

# ===== 细胞类型注释 =====
cat("\n===== 3. 细胞类型注释 =====\n")

# 基于marker基因表达和各cluster特征进行注释
# 根据分析结果，手动分配细胞类型
seurat_object$cell.types <- "Unknown"
seurat_object$cell.types[seurat_object$seurat_clusters %in% c(0, 3)] <- "Myeloid cells"
seurat_object$cell.types[seurat_object$seurat_clusters %in% c(1, 8)] <- "NK/T cells"
seurat_object$cell.types[seurat_object$seurat_clusters %in% c(2)] <- "Osteoblastic OS cells"
seurat_object$cell.types[seurat_object$seurat_clusters %in% c(4)] <- "CAFs"
seurat_object$cell.types[seurat_object$seurat_clusters %in% c(5)] <- "Plasma cells"
seurat_object$cell.types[seurat_object$seurat_clusters %in% c(6)] <- "Osteoclasts OCs"
seurat_object$cell.types[seurat_object$seurat_clusters %in% c(7)] <- "Endothelial cells"
seurat_object$cell.types[seurat_object$seurat_clusters %in% c(9)] <- "B cells"

# 统计各细胞类型数量
cat("\n细胞类型分布:\n")
print(table(seurat_object$cell.types))

# 可视化
color <- c("#e7cd79","#c1e6f3","#e69f84","#479d88","#415284",
           "#aa3538","#5891bf","#ab3282","#c5deba","#9fa3a8")
p <- DimPlot(seurat_object, group.by = "cell.types", reduction = "umap", label = TRUE,
             pt.size = 0.1, cols = color)
ggsave("graph/09_CellType.pdf", width = 8, height = 6, plot = p)
ggsave("graph/09_CellType.tiff", width = 8, height = 6, dpi = 300, plot = p)

# 保存结果
save(seurat_object, file = "seurat_object.rdata")

cat("\n===== 注释完成! =====\n")