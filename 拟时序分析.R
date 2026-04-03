# ===== 镜像设置（国内服务器使用）=====
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/BiocMirror/")

# 安装包（如需要）
# BiocManager::install("monocle", force = T)
# remotes::install_github("cole-trapnell-lab/monocle3")

library(monocle)
library(Seurat)
library(ggpubr)
library(monocle3)

#monocle 和 monocle3不能同时library

rm(list = ls())
gc()
setwd("E:/22_github_uploaded/Single _cell")

load("scRNA_NK&T.rdata")

# monocle3 拟时序分析 ----------------------------------------------------------
remotes::install_github("cole-trapnell-lab/monocle3")
library(monocle3)

# 添加行名和列名
rownames(scRNAsub@assays$RNA@layers$counts) = Features(scRNAsub)
colnames(scRNAsub@assays$RNA@layers$counts) = Cells(scRNAsub)
rownames(scRNAsub@assays$RNA@layers$data) = Features(scRNAsub)
colnames(scRNAsub@assays$RNA@layers$data) = Cells(scRNAsub)

# 获取表达矩阵
expr = GetAssayData(scRNAsub, assay = "RNA", layer = "counts")
# 获取表型信息
cell_metadata = scRNAsub@meta.data
# 获取基因名
gene_annotation = data.frame(gene_short_name = rownames(expr),
                             row.names = row.names(expr))
# 构建CDS对象
cds = new_cell_data_set(expr,
                        cell_metadata = cell_metadata,
                        gene_metadata = gene_annotation)

# NormalizeData+ScaleData+RunPCA
cds = preprocess_cds(cds,
                     method = "PCA",
                     num_dim = 50  #可以看图改成趋近于平衡的数目，也可以不改
                     )
plot_pc_variance_explained(cds)

# 去批次（可选）
cds = align_cds(cds, alignment_group = "sample")
# 降维
cds = reduce_dimension(cds,
                       preprocess_method = "PCA",  # "LSI","PCA",需要与preprocess_cds中一致
                       reduction_method = "UMAP")  # "UMAP","tSNE","PCA","LSI","Aligned"
# 可视化
p4 = plot_cells(cds, color_cells_by = "cell.types")
p4
#会自己生成一个umap结果，但不需要，我们使用之前seurat的umap结果

# 导入之前seurat降维的UMAP结果
cds.uamp = cds@int_colData$reducedDims$UMAP
sc.uamp = Embeddings(scRNAsub, reduction = "umap")
sc.uamp = sc.uamp[rownames(cds.uamp),]  #如果细胞数目有变化，使用这句命令将共有的筛选出来
cds@int_colData$reducedDims$UMAP = sc.uamp
# 可视化
p5 = plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cell.types")
p5
#是之前细胞亚群重聚类的结果

# 查看指定基因表达情况
genes = c("CD74")
plot_cells(cds, genes = genes, label_cell_groups = F, show_trajectory_graph = F)

# 分区
#对不同的分化轨迹分区
cds = cluster_cells(cds, cluster_method = "louvain") # 参数louvain可不用输入分辨率
# 可视化
plot_cells(cds, color_cells_by = "partition") 

# 构建细胞轨迹，对不同分区的细胞进行单独的轨迹分析 (运行速度较慢)
cds = learn_graph(cds)
# save(cds, file = "data/cds_learn_graph.rdata")
# load("D:/R/scRNA/Class3/data/cds_learn_graph.rdata")
# 可视化
plot_cells(cds,
           color_cells_by = "partition",
           label_groups_by_cluster = FALSE, 
           label_leaves = F,
           label_branch_points = F)

# 确定起始位点
# 自动计算确定起始位点
get_earliest_principal_node = function(cds, time_bin = "0") {
  cell_ids = which(colData(cds)[, "seurat_clusters"] == 0)
  closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes = igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))
# 可视化
p6 = plot_cells(cds,
                color_cells_by = "pseudotime",
                label_cell_groups = FALSE,
                label_leaves = FALSE,     # 是否为主图形中的每个叶节点绘制标签
                label_branch_points = FALSE,  # 是否为主图形中的每个分支点绘制标签
                graph_label_size = 1.5)  # 创建branch,root和leaf标签的大小
p6
ggsave("monocle3_plot_cells.pdf", width = 8, height = 6, plot = p6)
ggsave("monocle3_plot_cells.tiff", width = 8, height = 6, dpi = 300, plot = p6)

# 手动确定起始位点
cds = order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5)

# 某一基因表达情况可视化
p8 = monocle3::plot_genes_in_pseudotime(cds["CD74",],
                                        color_cells_by = "RNA_snn_res.0.2",
                                        min_expr = 0.5,
                                        ncol = 2)
p8
ggsave("monocle3_plot_genes1.pdf", width = 8, height = 6, plot = p8)
ggsave("monocle3_plot_genes1.tiff", width = 8, height = 6, dpi = 300, plot = p8)
p9 = monocle3::plot_genes_in_pseudotime(cds["CD74",],
                                        color_cells_by = "cell.types",
                                        min_expr = 0.5,
                                        ncol = 2)
p9
ggsave("monocle3_plot_genes2.pdf", width = 8, height = 6, plot = p9)
ggsave("monocle3_plot_genes2.tiff", width = 8, height = 6, dpi = 300, plot = p9)


