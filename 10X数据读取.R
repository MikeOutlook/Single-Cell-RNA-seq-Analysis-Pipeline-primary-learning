
# ===== 镜像设置（国内服务器使用）=====
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/BiocMirror/")

# 安装R包
if (!require("Seurat", quietly = TRUE)) install.packages("Seurat")
#单细胞分析包
# 加载R包
library(dplyr)
library(data.table)
library(Seurat)

rm(list = ls())
gc()
setwd("E:/22_github_uploaded/Single _cell")


# 设置10X数据存放路径
path_10X = "data"
# 获取上述路径下的所有样本文件列表
dir = list.files(path_10X)

##### 创建Seurat对象 #####
#先创建seurat数据再合并样本
# 创建一个空的列表来存储Seurat对象
seurat_list = list()
# 循环读取每个样本的10x数据并创建Seurat对象
for (i in 1:length(dir)) {
  print(i)
  data.path = file.path(path_10X, dir[i])  # 拼接文件路径
  seurat_data = Read10X(data.dir = data.path)
  # 创建Seurat对象
  seurat = CreateSeuratObject(counts = seurat_data,
                              project = dir[i],  # 指定项目名称为样本文件名
                              min.features = 200,# 每个细胞最少表达200个基因
                              min.cells = 3)     # 每个基因至少在多少个细胞中表达
  seurat_list = append(seurat_list, seurat)      # 将Seurat对象添加到列表中
}

# 合并Seurat对象，基因取并集，细胞加起来
seurat_object = merge(seurat_list[[1]], 
                      y = seurat_list[-1],
                      add.cell.ids = dir #添加细胞对应的样本信息，区分不同细胞位于哪个样本
                      )
#合并count数据
seurat_object = JoinLayers(seurat_object)
View(seurat_object@meta.data)
#meta.data中，nFeature代表细胞中表达的总基因数
#nCount代表一个基因中基因表达的总数
table(Idents(seurat_object))
#Idents指定细胞类型/细胞分组，例如T细胞，第一组数据等
seurat_object$sample = Idents(seurat_object)
#将样本信息储存到sample中


# Seurat对象保存 ---------------------------------------------------------------
save(seurat_object, file = "GSE261170_data.rdata")


