# ===== 镜像设置（国内服务器使用）=====
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/BiocMirror/")

library(Seurat)
library(ggplot2)
library(dplyr)

setwd("E:/22_github_uploaded/Single _cell")
load("seurat_object.rdata")

#添加行名列名
rownames(seurat_object@assays$RNA@layers$counts) = Features(seurat_object)
colnames(seurat_object@assays$RNA@layers$counts) = Cells(seurat_object)
rownames(seurat_object@assays$RNA@layers$data) = Features(seurat_object)
colnames(seurat_object@assays$RNA@layers$data) = Cells(seurat_object)

# 手动注释
# 方法一 -----------------------------------------------------------------------
#需要用到文章信息，大部分文章都会提供注释的分类以及marker
#根据提供的marker，找到这些marker在哪些cluster里高表达，综合考虑每个簇里哪组基因表达最高以及umap图来确定cluster的注释
markers = c("ALPL","RUNX2","IBSP",# Osteoblastic OS cells  2
            "LYZ","CD68",         # Myeloid cells  0 3
            "ACP5","CTSK",        # OCs  6
            "COL1A1","FAP","VIM", # CAFs  4
            "CD2","CD3D","CD3E","CD3G","GNLY","NKG7","KLRD1","KLRB1", # NK/T cells  1 8
            "EGFL7","PLVAP", # Endothelial cells  7
            "MS4A1","CD79A", # B cells  9
            "IGHG1","MZB1")  # Plasma cells  5
p1 = FeaturePlot(seurat_object, features = markers, ncol = 6)
p1
ggsave("graph/Markers1.tiff", width = 16, height = 12, dpi = 300, plot = p1)
p2 = DotPlot(seurat_object, features = markers) + theme_bw() + RotatedAxis()
p2
ggsave("graph/Markers2.tiff", width = 8, height = 6, dpi = 300, plot = p2)
p3 = VlnPlot(seurat_object, features = markers, stack = T, flip = T) + NoLegend()
p3
ggsave("graph/Markers3.tiff", width = 8, height = 12, dpi = 300, plot = p3)

# 根据聚类分布辅助判断
DimPlot(seurat_object, label = T, reduction = "umap")
# 添加细胞注释信息
seurat_object$cell.types = recode(seurat_object$seurat_clusters,
                                  "0" = "Myeloid cells",
                                  "1" = "NK/T cells",
                                  "2" = "Osteoblastic OS cells",
                                  "3" = "Myeloid cells",
                                  "4" = "CAFs",
                                  "5" = "Plasma cells",
                                  "6" = "OCs",
                                  "7" = "Endothelial cells",
                                  "8" = "NK/T cells",
                                  "9" = "B cells",
                                  "10" = "Unknown",
                                  "11" = "Unknown"
)
table(seurat_object$cell.types)
color = c("#e7cd79","#c1e6f3","#e69f84","#479d88","#415284",
          "#aa3538","#5891bf","#ab3282","#c5deba","#9fa3a8")
p1 = DimPlot(seurat_object, group.by = "cell.types", reduction = "umap", label = T, pt.size = 0.1, cols = color)
p1
ggsave("graph/09_CellType.pdf", width = 8, height = 6, plot = p1)
ggsave("graph/09_CellType.tiff", width = 8, height = 6, plot = p1, dpi = 300)
save(seurat_object, file = "seurat_object.rdata")


# 方法二 -----------------------------------------------------------------------
# 寻找marker，对每一个cluster进行运算，比较慢
allmarkers = FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#原理：0作为一组，其它簇作为一组进行差异分析，找到0中的marker；每一个簇进行以上操作找到每一个簇中的marker
#min.pct：基因在群体中的最小表达比例，低于该比例的基因会被过滤。
#logfc.threshold：差异倍数的绝对值
write.csv(allmarkers, "allmarkers.csv")
# allmarkers = read.csv("data/allmarkers.csv", row.names = 1, check.names = F)
#每一个簇选择前10个marker
Top10.coarse = allmarkers %>% 
  group_by(cluster) %>% 
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(Top10.coarse, "Top10.coarse.csv")
#cellmarker网站注释http://117.50.127.228/CellMarker/
#https://panglaodb.se/

seurat_object$cell.types = recode(seurat_object$seurat_clusters,
                                  "0" = "Myeloid cells",
                                  "1" = "Osteoblastic OS cells",
                                  "2" = "NK&T cells",
                                  "3" = "Myeloid cells",
                                  "4" = "Myeloid cells",
                                  "5" = "CAFs",
                                  "6" = "Plasma cells",
                                  "7" = "OCs",
                                  "8" = "Endothelial cells",
                                  "9" = "B cells",
                                  "10" = "Myeloid cells",
                                  "11" = "Plasma cells",
                                  "12" = "Osteoblastic OS cells")
table(seurat_object$cell.types)
color = c("#e7cd79","#c1e6f3","#e69f84","#479d88","#415284",
          "#aa3538","#5891bf","#ab3282","#c5deba","#9fa3a8")
p1 = DimPlot(seurat_object, group.by = "cell.types", reduction = "umap", label = T, pt.size = 0.1, cols = color)
p1
ggsave("graph/09_CellType.pdf", width = 8, height = 6, plot = p1)
ggsave("graph/09_CellType.tiff", width = 8, height = 6, plot = p1, dpi = 300)



# 自动注释 ---------------------------------------------------------------------
#不推荐
BiocManager::install("SingleR")
BiocManager::install("celldex")

library(SingleR)
library(celldex)

# singleR进行单细胞注释
# HumanPrimaryCellAtlasData(...) human 
# BlueprintEncodeData(...) human
# DatabaseImmuneCellExpressionData(...) human
# NovershternHematopoieticData(...) human
# MonacoImmuneData(...) human
# ImmGenData(...) mouse
# MouseRNAseqData(...) mouse

# 加载参考数据集
ref_Human_all = celldex::HumanPrimaryCellAtlasData()
save(human.ref, file = "SingleR_ref/ref_Human_all.RData")
load("SingleR_ref/ref_Human_all.RData")
load("SingleR_ref/ref_Hematopoietic.RData")

# 一个参考数据集进行细胞注释
cell.type = SingleR(test = as.matrix(seurat_object@assays[["RNA"]]@layers[["data"]]), # 输入表达矩阵
                    ref = ref_Human_all,                                              # 参考数据
                    labels = ref_Human_all$label.fine,                                # 标签列
                    clusters = seurat_object$seurat_clusters)                         # 细胞聚类信息
# 查看注释情况
table(cell.type$labels, levels(seurat_object$seurat_clusters))
cells = data.frame(id = levels(seurat_object$seurat_clusters), cell.type = cell.type$labels)
plotScoreHeatmap(cell.type)
# 加载注释信息
seurat_object@meta.data$labels = cells[match(seurat_object$seurat_clusters, cells$id), 2]

# 多个参考数据集进行注释
cell.type = SingleR(test = as.matrix(seurat_object@assays[["RNA"]]@layers[["data"]]), # 输入表达矩阵
                    ref = list(BP = ref_Hematopoietic, HPCA = ref_Human_all),         # 参考数据
                    # ref = list(a = ref1, b = ref2, c = ref3)
                    labels = list(ref_Hematopoietic$label.fine, ref_Human_all$label.fine),# 标签列
                    # labels = list(ref1$label.fine, ref2$label.fine, ref3$label.fine)
                    clusters = seurat_object$seurat_clusters)                         # 细胞聚类信息
# 查看注释情况
table(cell.type$labels, levels(seurat_object$seurat_clusters))
cells = data.frame(id = levels(seurat_object$seurat_clusters), cell.type = cell.type$labels)
plotScoreHeatmap(cell.type)
# 加载注释信息
seurat_object@meta.data$labels1 = cells[match(seurat_object$seurat_clusters, cells$id), 2]

# 比较注释信息
table(seurat_object$labels, seurat_object$labels1)
p1 = DimPlot(seurat_object, label = T)
p2 = DimPlot(seurat_object, group.by = "labels", label = T)
p1+p2



# 不同细胞类型占比----------------------------
metadata=seurat_object[[]]
# metadata$sample=factor(metadata$sample,levels = c("sham","CLP1","CLP2"))

color = c("#e7cd79","#c1e6f3","#e69f84","#479d88","#415284",
          "#aa3538","#5891bf","#ab3282","#c5deba","#9fa3a8","#9370DB")
plot_group<-ggplot(metadata,aes(x=sample,fill=cell.types))+
  geom_bar(position="fill")+
  scale_fill_manual(values=color) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")
plot_group
ggsave("graph/Celltypes_pro.pdf", width = 8, height = 6, plot = plot_group)
ggsave("graph/Celltypes_pro.tiff", width = 8, height = 6, plot = plot_group, dpi = 300)












