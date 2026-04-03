
# 安装R包
install.packages("Seurat")
devtools::install_github(repo = "hhoeflin/hdf5r")
# 加载R包
library(dplyr)
library(data.table)
library(Seurat)
library(hdf5r)

rm(list = ls())
gc()
setwd("/media/desk16/brs7/WZH/class2/h5")


# 设置h5数据存放路径
path_h5 = "/media/desk16/brs7/WZH/class2/h5/GSE200874"
# 获取上述路径下的所有样本文件列表
dir = list.files(path_h5, pattern = "\\.h5$")
#设置样本ID，从文件名中提取或者自定义
sample = sapply(strsplit(dir, "_"), function(x) {x[[1]]})
sample = c("sample1","sample2","sample3","sample4")

##### 创建Seurat对象 #####
# 创建一个空的列表来存储Seurat对象
seurat_list = list()
# 循环读取每个样本的10x数据并创建Seurat对象
for (i in 1:length(dir)) {
  print(i)
  data.path = file.path(path_h5, dir[i])  # 拼接文件路径
  seurat_data = Read10X_h5(filename = data.path)
  # 创建Seurat对象
  seurat = CreateSeuratObject(counts = seurat_data,
                              project = sample[i],# 指定项目名称为样本文件名
                              min.features = 200, # 每个细胞最少表达200个基因
                              min.cells = 3)      # 每个基因至少在多少个细胞中表达
  seurat_list = append(seurat_list, seurat)       # 将Seurat对象添加到列表中
}

# 合并Seurat对象
seurat_object = merge(seurat_list[[1]], 
                      y = seurat_list[-1],
                      add.cell.ids = sample)
seurat_object = JoinLayers(seurat_object)
table(Idents(seurat_object))
seurat_object$sample = Idents(seurat_object)



# Seurat对象保存 ---------------------------------------------------------------
save(seurat_object, file = "/media/desk16/brs7/WZH/class2/h5/h5_data.rdata")


