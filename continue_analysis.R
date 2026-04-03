# ===== 重新创建Seurat对象并继续分析 =====
options(repos = c(CRAN = "https://cloud.r-project.org/"))

library(dplyr)
library(data.table)
library(ggplot2)
library(Seurat)
library(tidyverse)
library(clustree)
library(patchwork)
library(harmony)

rm(list = ls())
gc()
options(future.globals.maxSize = 10*1024^3)
setwd("/root/autodl-tmp/Single _cell")

cat("===== 1. 读取10X数据并创建Seurat对象 =====\n")
path_10X <- "data"
dir <- list.files(path_10X)
cat("找到", length(dir), "个样本\n")

seurat_list <- list()
for (i in 1:length(dir)) {
  cat("处理样本", i, "/", length(dir), ":", dir[i], "\n")
  data.path <- file.path(path_10X, dir[i])
  seurat_data <- Read10X(data.dir = data.path)
  seurat <- CreateSeuratObject(counts = seurat_data,
                              project = dir[i],
                              min.features = 200,
                              min.cells = 3)
  seurat_list <- append(seurat_list, seurat)
}

seurat_object <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = dir)
# Seurat 4.x 不需要 JoinLayers

cat("总细胞数:", ncol(seurat_object), "\n")
cat("总基因数:", nrow(seurat_object), "\n")

# 添加分组信息
seurat_object$sample <- Idents(seurat_object)
seurat_object$condition <- as.vector(seurat_object$sample)
seurat_object$condition[which(seurat_object$condition=="GSM8136885")] <- "EarlyAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136886")] <- "PermanentAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136887")] <- "PermanentAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136888")] <- "EarlyAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136889")] <- "EarlyAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136890")] <- "EarlyAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136897")] <- "PermanentAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136898")] <- "PermanentAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136899")] <- "PermanentAF"

cat("\n===== 2. 质控 =====\n")
# 线粒体等比例
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
seurat_object[["percent.hb"]] <- PercentageFeatureSet(seurat_object, pattern = "^HB[^(P)]")
seurat_object[["percent.ribo"]] <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]")

# 质控过滤
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 10)
cat("质控后细胞数:", ncol(seurat_object), "\n")

cat("\n===== 3. 归一化 =====\n")
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

cat("\n===== 4. 细胞周期评分 =====\n")
seurat_object <- CellCycleScoring(seurat_object,
                               s.features = cc.genes$s.genes,
                               g2m.features = cc.genes$g2m.genes)

cat("\n===== 5. 缩放和批次校正 =====\n")
# Seurat 4.x 直接缩放，不需要 split
seurat_object <- ScaleData(seurat_object, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))

cat("\n===== 6. PCA =====\n")
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

cat("\n===== 7. Harmony去批次 =====\n")
seurat_object <- RunHarmony(seurat_object, group.by.vars = "sample", plot_convergence = TRUE)
# Seurat 4.x 不需要 JoinLayers

cat("\n===== 8. 聚类 =====\n")
seurat_object <- FindNeighbors(seurat_object, reduction = "harmony", dims = 1:20)
seurat_object <- FindClusters(seurat_object, resolution = seq(from = 0.1, to = 1.0, by = 0.1))

cat("\n===== 9. 分辨率0.2 =====\n")
seurat_object <- FindClusters(seurat_object, resolution = 0.2)
Idents(seurat_object) <- "RNA_snn_res.0.2"

cat("\n===== 10. UMAP和TSNE =====\n")
seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = 1:20)
seurat_object <- RunTSNE(seurat_object, reduction = "harmony", dims = 1:20)

cat("\n===== 11. 保存结果 =====\n")
save(seurat_object, file = "seurat_object.rdata")
cat("分析完成! 已保存为 seurat_object.rdata\n")
cat("最终细胞数:", ncol(seurat_object), "\n")
cat("cluster数:", length(levels(seurat_object@active.ident)), "\n")