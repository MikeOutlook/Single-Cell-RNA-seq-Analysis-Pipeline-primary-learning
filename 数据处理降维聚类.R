
# ===== 镜像设置（国内服务器使用）=====
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/BiocMirror/")

# R包已在之前安装，这里跳过安装步骤
# install.packages(c("ggplot2", "Seurat", "tidyverse", "clustree", "patchwork"))
# BiocManager::install("decontX")
# devtools::install_github("immunogenomics/harmony")

library(dplyr)
library(data.table)
library(ggplot2)
library(Seurat)
library(tidyverse)
library(clustree)
# library(decontX)  # 跳过decontX
library(patchwork)
library(harmony)

rm(list = ls())
gc()
options(future.globals.maxSize = 10*1024^3)
setwd("E:/22_github_uploaded/Single _cell")

# 读取seurat对象
load("GSE261170_data.rdata")
View(seurat_object@meta.data)
table(seurat_object$sample)

# 给counts数据添加行名和列名
rownames(seurat_object@assays$RNA@layers$counts) = Features(seurat_object)
colnames(seurat_object@assays$RNA@layers$counts) = Cells(seurat_object)

# 添加分组信息
# GSE261170 样本分组：
# Early AF: GSM8136885, GSM8136888, GSM8136889, GSM8136890
# Permanent AF: GSM8136886, GSM8136887, GSM8136897, GSM8136898, GSM8136899
seurat_object$condition = as.vector(seurat_object$sample)
seurat_object$condition[which(seurat_object$condition=="GSM8136885")] = "EarlyAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136886")] = "PermanentAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136887")] = "PermanentAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136888")] = "EarlyAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136889")] = "EarlyAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136890")] = "EarlyAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136897")] = "PermanentAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136898")] = "PermanentAF"
seurat_object$condition[which(seurat_object$condition=="GSM8136899")] = "PermanentAF"
table(seurat_object$sample)
table(seurat_object$condition)

# 创建线粒体,红细胞,核糖体基因表达比例
View(seurat_object[[]]) #提取meta信息
seurat_object[["percent.mt"]] = PercentageFeatureSet(seurat_object, pattern = "^MT-")  #MT开头的都是线粒体基因，计算线粒体基因占比
seurat_object[["percent.hb"]] = PercentageFeatureSet(seurat_object, pattern = "^HB[^(P)]")  #红细胞
seurat_object[["percent.ribo"]] = PercentageFeatureSet(seurat_object, pattern = "^RP[SL]")  #核糖体基因表达比例
#若是小鼠的数据，MT/HB/RP改为小写
# 可视化,便于后续数据质控
p1 = VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 展示红细胞和核糖体的图，在features中加percent.hb,percent.ribo
p1
#质控之前的样本信息
ggsave("graph/01_QCbefore.pdf", width = 15, height = 6, plot = p1)
ggsave("graph/01_QCbefore.tiff", width = 15, height = 6, dpi = 300, plot = p1)

# 相关性分析，判断数据质量是否可以
plot1 = FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")  #该图相关性越趋近于0越好，越大说明线粒体基因表达越多，细胞濒临死亡
plot2 = FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  #相关性系数越趋近于1越好
p2 = plot1 + plot2
p2
ggsave("graph/02_FeatureScatter.pdf", width = 12, height = 6, plot = p2)
ggsave("graph/02_FeatureScatter.tiff", width = 12, height = 6, dpi = 300, plot = p2)

# 数据质控 针对线粒体相关基因表达比例以及nFeature值进行质控
#筛选范围参考原始文章
#大部分文章不筛选红细胞和核糖体
seurat_object = subset(seurat_object, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 10)
table(seurat_object$orig.ident)
# 质控后可视化
p3 = VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
p3
ggsave("graph/03_QCafter.pdf", width = 15, height = 6, plot = p3)
ggsave("graph/03_QCafter.tiff", width = 15, height = 6, dpi = 300, plot = p3)

# 数据全局缩放归一化
seurat_object = NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000) #标准化，数据在count数据下新生成一个data数据
seurat_object@assays$RNA@layers$counts[1:15,1:15] # 稀疏矩阵 点代表该基因在对应的细胞中表达值为0，为了让内存更小，0显示为点
seurat_object@assays$RNA@layers$data[1:15,1:15]

# 寻找高变基因，便于后续聚类
seurat_object = FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
# 查看前10个高变基因选择情况
top10 = head(VariableFeatures(seurat_object), 10)
plot1 = VariableFeaturePlot(seurat_object)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0) 
combined_plot = plot1 + plot2 + plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "bottom")
print(combined_plot)
#红色的点为高变基因，黑色的为非高变基因

# 估算细胞周期效应
# 上一步找到的高变基因中常常会包含一些细胞周期相关基因
# 它们会导致细胞聚类发生一定的偏移,即相同类型的细胞在聚类时会因为细胞周期的不同而分开
# 细胞周期分为s期，G1期，G2期和M期，这里把G2期和M期放在一起了
# CaseMatch(c(cc.genes$s.genes, cc.genes$g2m.genes), VariableFeatures(seurat_object)) #与细胞周期有关的基因
g2m_genes = CaseMatch(search = cc.genes$g2m.genes, match = rownames(seurat_object)) #提取g2m期基因
s_genes = CaseMatch(search = cc.genes$s.genes, match = rownames(seurat_object))  #提取s期基因
seurat_object = CellCycleScoring(seurat_object,s.features = cc.genes$s.genes,
                                 g2m.features = cc.genes$g2m.genes)  #计算评分
#细胞周期信息存储在meta数据中
# 查看细胞周期基因表达情况
p4 = RidgePlot(seurat_object, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)  #这四个基因是前人研究中对细胞周期更有代表性的基因
p4
# 峰值在0表示不受细胞周期影响，若峰值不在0表示收到细胞周期的影响，需要矫正
ggsave("graph/04_RidgePlot_CellCycleScore.pdf", width = 8, height = 6, plot = p4)
ggsave("graph/04_RidgePlot_CellCycleScore.tiff", width = 8, height = 6, dpi = 300, plot = p4)


VlnPlot(seurat_object,features = "TOP2A")


# 去除游离RNA污染 - 跳过decontX
# decontX_results = decontX::decontX(seurat_object@assays$RNA@layers$counts)
# seurat_object$Contamination = decontX_results$contamination

# 校正线粒体影响(和细胞周期效应影响)
#根据sample分割出不同组数据的sample
seurat_object[["RNA"]] = split(seurat_object[["RNA"]], f = seurat_object$sample)
#归一化处理，同时矫正线粒体对数据的影响
seurat_object = ScaleData(seurat_object, vars.to.regress = c("percent.mt"))
#矫正线粒体和细胞周期对数据的影响
seurat_object = ScaleData(seurat_object, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))

# PCA降维
seurat_object = RunPCA(seurat_object, features = VariableFeatures(object = seurat_object)) #基于之前选择的2000个高变基因
DimPlot(object = seurat_object, reduction = "pca", pt.size = 0.1, group.by = "sample")
#若每个样本的点都混合在一块区域，说明该数据没有批次效应
#该数据有批次效应

# 去批次:去除样品间的差异
#一般两个都要运行，选择更好的一个
# 方法一 使用CCA去批次
# 运行很久
seurat_object = IntegrateLayers(object = seurat_object,
                                method = CCAIntegration,
                                orig.reduction = "pca", #基于PCA去批次
                                new.reduction = "cca",  #生成cca用于存数据
                                verbose = FALSE  #显示进度条
                                )
p5 = DimPlot(object = seurat_object, reduction = "cca", pt.size = 0.1, group.by = "sample")
p5
ggsave("graph/05_CCA.pdf", width = 8, height = 6, plot = p5)
ggsave("graph/05_CCA.tiff", width = 8, height = 6, dpi = 300, plot = p5)

# 方法二 使用Harmony去批次
seurat_object = RunHarmony(seurat_object, group.by.vars = "sample", plot_convergence = TRUE)
p6 = DimPlot(object = seurat_object, reduction = "harmony", pt.size = 0.1, group.by = "sample")
p6
ggsave("graph/05_Harmony.pdf", width = 8, height = 6, plot = p6)
ggsave("graph/05_Harmony.tiff", width = 8, height = 6, dpi = 300, plot = p6)

# 合并矩阵
seurat_object = JoinLayers(seurat_object)

# 降维可视化
# 查看最合适的主成分个数，pc数选择尽量小的，一般20-30之间
p7 = ElbowPlot(seurat_object, ndims = 50, reduction = "pca")+theme_bw()
p7
ggsave("graph/06_Elbowplot.pdf", width = 8, height = 6, plot = p7)
ggsave("graph/06_Elbowplot.tiff", width = 8, height = 6, dpi = 300, plot = p7)

# 进行谱聚类,基于SNN和模块化优化的聚类算法来识别细胞簇
seurat_object = FindNeighbors(seurat_object, reduction = "harmony", dims = 1:20) #dims参数写pc数
# 调整分辨率,当我们样本细胞数较大时resolution要高一些
#分辨率越高，cluster越多
#从0.1到1分别计算一个结果，文章中若有提到，直接设置resolution参数即可
seurat_object = FindClusters(seurat_object, resolution = seq(from = 0.1, to = 1.0, by = 0.1))
p8 = clustree(seurat_object)
p8
#选择分辨率，浅色箭头越少越好
ggsave("graph/07_clustree.pdf", width = 12, height = 10, plot = p8)
ggsave("graph/07_clustree.tiff", width = 12, height = 10, dpi = 300, plot = p8)

seurat_object = FindClusters(seurat_object, resolution = 0.2)
Idents(seurat_object) = "RNA_snn_res.0.2"
# 聚类可视化
seurat_object = RunUMAP(seurat_object, reduction = "harmony", dims = 1:20)
seurat_object = RunTSNE(seurat_object, reduction = "harmony", dims = 1:20)
p9 = DimPlot(seurat_object, reduction = "umap", label = T)
p10 = DimPlot(seurat_object, reduction = "tsne", label = T)
p9+p10
ggsave("graph/08_UMAP.pdf", width = 8, height = 6, plot = p9)
ggsave("graph/08_UMAP.tiff", width = 8, height = 6, dpi = 300, plot = p9)
ggsave("graph/08_TSNE.pdf", width = 8, height = 6, plot = p10)
ggsave("graph/08_TSNE.tiff", width = 8, height = 6, dpi = 300, plot = p10)

# 保存Seurat对象
save(seurat_object, file = "seurat_object.rdata")


