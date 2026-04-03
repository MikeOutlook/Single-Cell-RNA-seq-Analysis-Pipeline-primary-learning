# ===== 镜像设置（国内服务器使用）=====
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/BiocMirror/")

# install.packages(c("dplyr", "data.table", "ggplot2", "Seurat", "tidyverse", "clustree", "patchwork", "harmony"))

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
setwd("E:/22_github_uploaded/Single _cell")

# 细胞亚类再聚类 ---------------------------------------------------------------
#数据已经进行过去批次，没有再次去批次
load("seurat_object.rdata")
table(seurat_object$cell.types)
scRNAsub = subset(seurat_object, cell.types == "NK/T cells")
rownames(scRNAsub@assays$RNA@layers$counts) = Features(scRNAsub)
colnames(scRNAsub@assays$RNA@layers$counts) = Cells(scRNAsub)
rownames(scRNAsub@assays$RNA@layers$data) = Features(scRNAsub)
colnames(scRNAsub@assays$RNA@layers$data) = Cells(scRNAsub)
# 标准化
scRNAsub = NormalizeData(scRNAsub, normalization.method = "LogNormalize", scale.factor = 10000)  #标准化处理
# 寻找高变基因
scRNAsub = FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
# 归一化
scRNAsub[["RNA"]] = split(scRNAsub[["RNA"]], f = scRNAsub$sample)
scRNAsub = ScaleData(scRNAsub, vars.to.regress = c("percent.mt"))
# PCA
scRNAsub = RunPCA(scRNAsub, features = VariableFeatures(object = scRNAsub))
DimPlot(object = scRNAsub, reduction = "pca", pt.size = 0.1, group.by = "sample")
# # 使用CCA去批次
# scRNAsub = IntegrateLayers(object = scRNAsub,
#                            method = CCAIntegration,
#                            orig.reduction = "pca",
#                            new.reduction = "cca",
#                            verbose = FALSE)
# DimPlot(object = scRNAsub, reduction = "cca", pt.size = 0.1, group.by = "sample")
# 使用Harmony去批次
scRNAsub = RunHarmony(scRNAsub, group.by.vars = "sample", plot_convergence = TRUE)
DimPlot(object = scRNAsub, reduction = "harmony", pt.size = 0.1, group.by = "sample")
# 合并矩阵
scRNAsub = JoinLayers(scRNAsub)

ElbowPlot(scRNAsub, ndims = 20, reduction="pca")
pc.num = 1:10
scRNAsub = FindNeighbors(scRNAsub, dims = pc.num, reduction = "harmony")
scRNAsub = FindClusters(scRNAsub, resolution = seq(from = 0.1, to = 1.0, by = 0.1))
clustree(scRNAsub)
scRNAsub$seurat_clusters = scRNAsub$RNA_snn_res.0.2
Idents(scRNAsub) = "RNA_snn_res.0.2"
scRNAsub = RunUMAP(scRNAsub, dims = pc.num,reduction = "harmony")
scRNAsub = RunTSNE(scRNAsub, reduction = "harmony", dims = pc.num)
p1=DimPlot(scRNAsub, reduction = "umap", label = T)
p2=DimPlot(scRNAsub, reduction = "tsne", label = T)
ggsave("immune_harmony.pdf", width = 8, height = 6, plot = p1)
ggsave("immune_harmony.tiff", width = 8, height = 6, dpi = 300, plot = p1)

#细胞注释
#检索T细胞亚群的marker基因
markers = c("HPGDS","KIT","FCER1A", # Mast cells  7
            "LAG3","TIGIT","CTLA4","CD2", # Tregs  3
            "NKG7","KLRD1","GZMB","GNLY", # NK cells  2
            "CD4","CD3E","CD3D","CD3G", # CD4-CD8- T cells   1
            "CD8A","CD8B") # CD8+ T cells  0  4  5  6
p2 = DotPlot(scRNAsub, features = markers) + RotatedAxis()
p2

# 根据聚类分布辅助判断
DimPlot(scRNAsub, label = T, reduction = "umap")
# 添加细胞注释信息
scRNAsub$cell.types = recode(scRNAsub$seurat_clusters,
                             "0" = "CD4-CD8- T cells",
                             "1" = "CD8+ T cells",
                             "2" = "CD8+ T cells",
                             "3" = "NK cells",
                             "4" = "Tregs",
                             "5" = "CD4-CD8- T cells",
                             "6" = "CD8+ T cells",
                             "7" = "Mast cells")
table(scRNAsub$cell.types)
color = c("#e7cd79","#c1e6f3","#e69f84","#479d88","#415284",
          "#aa3538","#5891bf","#ab3282","#c5deba","#9fa3a8")
p1 = DimPlot(scRNAsub, group.by = "cell.types", reduction = "umap", label = T, pt.size = 0.1, cols = color)
p1
ggsave("01_CellType.pdf", width = 8, height = 6, plot = p1)
ggsave("01_CellType.tiff", width = 8, height = 6, plot = p1, dpi = 300)
save(scRNAsub, file = "scRNA_NK&T.rdata")
