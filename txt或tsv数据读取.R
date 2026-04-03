
# 安装R包
install.packages("Seurat")
# 加载R包
library(dplyr)
library(data.table)
library(Seurat)

rm(list = ls())
gc()
setwd("/media/desk16/brs7/WZH/class2/matrix/")

# 有多个单细胞数据时 -----------------------------------------------------------

# 设置数据存放路径
path_txt = "/media/desk16/brs7/WZH/class2/matrix/GSE158035"
# 获取上述路径下的所有样本文件列表
dir = list.files(path_txt)
sample = sapply(strsplit(dir,"[.]"), function(x) {x[[1]]})
sample = sapply(strsplit(sample,"_"), function(x) {x[[2]]})

##### 创建Seurat对象 #####
# 创建一个空的列表来存储Seurat对象
seurat_list = list()
# 循环读取每个样本的10x数据并创建Seurat对象
for (i in 1:length(dir)) {
  print(i)
  data.path = file.path(path_txt, dir[i])  # 拼接文件路径
  seurat_data = read.table(gzfile(data.path), header = T, row.names = 1)
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
save(seurat_object, file = "/media/desk16/brs7/WZH/class2/matrix/txt_data.rdata")


