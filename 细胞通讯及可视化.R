#细胞通讯必须保证cell.types为factor格式，同时保证没有细胞个数为0的细胞类型

# ===== 镜像设置（国内服务器使用）=====
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/BiocMirror/")

# 安装包（如需要）
# install.packages(c("NMF", "tidydr"))
# devtools::install_github("jokergoo/circlize")
# devtools::install_github("sqjin/CellChat")

library(Seurat)
library(dplyr)
library(magrittr)
library(CellChat)
library(patchwork)
library(tidydr)

rm(list = ls())
gc()
setwd("E:/22_github_uploaded/Single _cell")

load("seurat_object.rdata")
metadata=seurat_object@meta.data

#疾病和正常细胞必须分开分析
table(metadata$condition)
# 疾病组cellchat分析 -----------------------------------------------------------
# 提取Case组的Seurat对象
case.object = subset(seurat_object, condition == "EarlyAF")
# 提取data数据
case.data = GetAssayData(case.object, assay = "RNA", layer = "data")
# 提取部分meta.data,主要是cell.types
case.meta = case.object@meta.data[,c("cell.types", "condition")]

# 构建cellchat对象,需要两个数据
case.cellchat = createCellChat(object = case.data)
case.cellchat = addMeta(case.cellchat, meta = case.meta)
# 将cell.types设置为默认的cell标识
case.cellchat = setIdent(case.cellchat, ident.use = "cell.types")
# 查看细胞类型
levels(case.cellchat@idents)

# cellchat提供的人类的配受体数据库
showDatabaseCategory(CellChatDB.human) #数据库基本信息
case.cellchat@DB = CellChatDB.human #小鼠导入CellChatDB.mouse
#如果有指定想要研究的配受体
# CellChatDB.use = subsetDB(CellChatDB.human, search = "Secreted Signaling", key="annotation")
# case.cellchat@DB = CellChatDB.use

# Secreted Signaling  专注于细胞间通讯，信号分子通过分泌途径从一个细胞释放到细胞外，然后影响附近的细胞，从而调节各种生物过程
# ECM-Receptor  主要关注细胞和基质的互动，信息可能更侧重于细胞在微环境中的位置和结构作用
# Cell-Cell Contact  主要关注细胞膜表面的直接接触和相互作用
# Heterodimers  专注于异二聚体的组合和功能，对受体和信号分子的具体组合提供信息
# KEGG  提供了广泛的生物通路图，包括细胞通讯，代谢通路，信号转导等，它是一个综合性数据库，覆盖了多种生物过程和通路
# Literature  提供的是研究成果和实验数据

# 表达数据预处理
case.cellchat = subsetData(case.cellchat) # 提取signaling gene,减少计算成本，提取配受体数据库中提及的基因
future::plan("multisession", workers = 12) # 设置线程数
case.cellchat = identifyOverExpressedGenes(case.cellchat) # 寻找高表达基因  结果输出在var.features里
case.cellchat = identifyOverExpressedInteractions(case.cellchat) # 寻找高表达基因的相互作用(通路)
case.cellchat = projectData(case.cellchat, PPI.human) # 将基因表达数据投射到PPI上  将基因层面的通路（相互作用）投射到蛋白层面，非必须

# cellchat分析
options(future.globals.maxSize = 2000 * 1024^2)
case.cellchat = computeCommunProb(case.cellchat, raw.use = T) # 使用表达值推测细胞互作的概率
# load("cellchat_Case.rdata")
#raw.use=T 基于data数据中的基因推测细胞互作概率
#将基因表达数据投射到 PPI 上（可选：运行时，用户应在函数 computeCommunProb() 中设置 raw.use = FALSE 以使用投射数据）
#raw.use=T对应基因数据，raw.use=F对应蛋白数据（PPI）
#raw.use=T 结果少一些但更容易与文献比较
case.cellchat = filterCommunication(case.cellchat, min.cells = 10) # 过滤掉小于10个细胞的胞间通讯网络
case.cellchat = computeCommunProbPathway(case.cellchat) # 计算每个信号通路相关的所有配体-受体相互作用的通信结果
case.cellchat = aggregateNet(case.cellchat) # 计算整合的细胞类型之间通信结果(配受体对数目和强度)
case.cellchat = netAnalysis_computeCentrality(case.cellchat, slot.name = "netP") # 计算网络中心性权重。可视化，判断哪个细胞类型是发射信号的一方，哪个细胞类型是接收信号的一方

# 数据查看,保存
Case.net = subsetCommunication(case.cellchat)
df.net <- subsetCommunication(case.cellchat, sources.use = c("OCs"), targets.use = c("CAFs")) #提取两个细胞之间的结果
df.net <- subsetCommunication(case.cellchat, signaling = c("WNT", "TGFb"))#提取部分信号通路结果

#每种细胞两两之间涉及到的配受体对
#source:发出信号的细胞类型
#target：接收信号的细胞类型
#ligand：配体
#receptor：受体
#prob：权重
write.csv(Case.net, "Case_net_inter_raw.useT.csv")
save(case.cellchat, file = "cellchat_Case.rdata")



# 正常组cellchat分析 -----------------------------------------------------------
# 提取Control组的Seurat对象
control.object = subset(seurat_object, condition == "PermanentAF")
# 提取data数据
control.data = GetAssayData(control.object, assay = "RNA", layer = "data")
# 提取部分meta.data,主要是cell.types
control.meta = control.object@meta.data[,c("cell.types", "condition")]

# 构建cellchat对象,需要两个数据
control.cellchat = createCellChat(object = control.data)
control.cellchat = addMeta(control.cellchat, meta = control.meta)
# 将cell.types设置为默认的cell标识
control.cellchat = setIdent(control.cellchat, ident.use = "cell.types")
# 查看细胞类型
levels(control.cellchat@idents) 
# 每种细胞类型对应的细胞数目
groupsize = as.numeric(table(control.cellchat@idents))
groupsize

# cellchat提供的人类的配受体数据库
control.cellchat@DB = CellChatDB.human #小鼠导入CellChatDB.mouse
# 表达数据预处理
control.cellchat = subsetData(control.cellchat) # 提取signaling gene,减少计算成本
future::plan("multisession", workers = 4) # 设置线程数
control.cellchat = identifyOverExpressedGenes(control.cellchat) # 寻找高表达基因
control.cellchat = identifyOverExpressedInteractions(control.cellchat) # 寻找高表达基因的相互作用(通路)
control.cellchat = projectData(control.cellchat, PPI.human) # 将基因表达数据投射到PPI上

# cellchat分析
control.cellchat = computeCommunProb(control.cellchat, raw.use = T) # 使用表达值推测细胞互作的概率
control.cellchat = filterCommunication(control.cellchat, min.cells = 10) # 过滤掉小于10个细胞的胞间通讯网络
control.cellchat = computeCommunProbPathway(control.cellchat) # 计算每个信号通路相关的所有配体-受体相互作用的通信结果
control.cellchat = aggregateNet(control.cellchat) # 计算整合的细胞类型之间通信结果(配受体对数目和强度)
control.cellchat = netAnalysis_computeCentrality(control.cellchat, slot.name = "netP") # 计算网络中心性权重

# 数据查看,保存
control.net = subsetCommunication(control.cellchat)
write.csv(control.net, "Control_net_inter_raw.useT.csv")
save(control.cellchat, file = "cellchat_Control.rdata")


# 结果可视化 -------------------------------------------------------------------
load("cellchat_Case.rdata")
load("cellchat_Control.rdata")
# 以疾病组为例
# 互作网络
cellchat = case.cellchat
groupsize = as.numeric(table(cellchat@idents))
#提取每个细胞类型包含的细胞数目


#p1:线条粗细代表互作的两种细胞类型的配受体数目
#p1:线条粗细代表互作的两种细胞类型的互作强度
pdf("graph/netVisual_circle.pdf", width = 12, height = 8)
par(mfrow = c(1,2), xpd = T)
p1 = netVisual_circle(cellchat@net$count,
                      vertex.weight = groupsize,
                      weight.scale = T,
                      label.edge = F,
                      title.name = "Number of interaction")
p2 = netVisual_circle(cellchat@net$weight,
                      vertex.weight = groupsize,
                      weight.scale = T,
                      label.edge = F,
                      title.name = "Interaction weights/strength")
dev.off()

# 检查每种细胞发出的信号
pdf("cellchat_every.pdf", width = 15, height = 12)
mat = cellchat@net$weight   # 提取每个细胞的互作强度
par(mfrow = c(3,3), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 = matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] = mat[i, ]
  netVisual_circle(mat2,
                   vertex.weight = groupsize,
                   weight.scale = T,
                   edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])
}
dev.off()
# 自定义可视化
levels(case.cellchat@idents)
par(mfrow = c(1,2), mar = c(1,1,1,1))
mat2 = matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#数字代表不同细胞类型，具体代表哪一种可以使用levels(case.cellchat@idents)查看
mat2[1, ] = mat[1, ]
netVisual_circle(mat2,
                 vertex.weight = groupsize,
                 weight.scale = T,
                 edge.weight.max = max(mat),
                 title.name = rownames(mat)[1])


# 选择感兴趣的通路可视化
# 层级图（Hierarchy plot）
cellchat@netP$pathways  #查看通路
pathways = c("MIF")
levels(cellchat@idents)
vertex.receiver = seq(1, 4)
#左边展示4个细胞类型
#线条粗细代表互作强度
pdf("graph/cellchat_Hierarchy.pdf", width = 15, height = 10)
netVisual_aggregate(cellchat,
                    signaling = pathways,
                    layout = "hierarchy",
                    vertex.receiver = vertex.receiver)
dev.off()
# 网络图（Circle plot）
pdf("graph/cellchat_circle.pdf", width = 6, height = 6)
netVisual_aggregate(cellchat,
                    signaling = pathways,
                    layout = "circle")
dev.off()
# 和弦图（Chord diagram）
pdf("cellchat_chord.pdf", width = 6, height = 6)
netVisual_aggregate(cellchat,
                    signaling = pathways,
                    layout = "chord")
dev.off()
# 热图（Heatmap）
pdf("graph/cellchat_Reds.pdf", width = 5, height = 4)
par(mfrow = c(1,1))
netVisual_heatmap(cellchat,
                  signaling = pathways,
                  color.heatmap = "Reds")
dev.off()


# 计算各个配受体对对指定信号通路的贡献
pdf("graph/contribution.pdf", width = 6, height = 3)
netAnalysis_contribution(cellchat, signaling = pathways)
dev.off()

# 提取对指定通路有贡献的所有配受体
all_LR = extractEnrichedLR(cellchat, signaling = pathways, geneLR.return = FALSE)
# 提取对这个通路贡献最大的配体受体对来展示（也可以选择其他的配体受体对）
LR = all_LR[1,]
# 层级图（Hierarchy plot）
vertex.receiver = seq(1,4)
pdf("graph/cellchat_Hierarchy_LR.pdf", width = 10, height = 8)
netVisual_individual(cellchat,
                     #color.use = cell_type_cols,
                     signaling = pathways,
                     pairLR.use = LR,
                     layout = "hierarchy",
                     vertex.receiver = vertex.receiver)
dev.off()
# 网络图（Circle plot）
pdf("graph/cellchat_circle_LR.pdf", width = 6, height = 6)
netVisual_individual(cellchat,
                     signaling = pathways,
                     pairLR.use = LR,
                     layout = "circle")
dev.off()
# 和弦图（Chord diagram）
pdf("cellchat_chord_LR.pdf", width = 6, height = 6)
netVisual_individual(cellchat,
                     signaling = pathways,
                     pairLR.use = LR,
                     layout = "chord")
dev.off()

# 多个配体-受体介导的细胞互作关系可视化----
table(cellchat@meta$cell.types)
# 气泡图，需要指定受体细胞和配体细胞
netVisual_bubble(cellchat, sources.use = c(1,3), targets.use = c(5:6), remove.isolate = FALSE)
# 气泡图（指定信号通路或配体-受体）
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:6), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
# 参与某条信号通路（如TGFb）的所有基因在细胞群中的表达情况展示（小提琴图）
plotGeneExpression(cellchat, signaling = "SPP1")


# 指定配受体细胞和通路(气泡图)
levels(cellchat@idents)
cellchat@netP$pathways  #查看通路
pathways = c("SPP1")
pdf("graph/netVisual_bubble.pdf", width = 10, height = 6)
netVisual_bubble(cellchat,
                 sources.use = c("OCs"),  #研究特定细胞的通路
                 targets.use = c("Myeloid cells","Osteoblastic OS cells","NK/T cells","CAFs"),
                 signaling = pathways,
                 remove.isolate = FALSE,
                 font.size = 14)
dev.off()


# 通过计算每个细胞群的网络中心性指标,识别每类细胞在信号通路中的角色/作用C（发送者、接收者、调解者和影响者）
netAnalysis_signalingRole_network(cellchat, signaling = pathways, width = 8, height = 2.5, font.size = 10)
# 识别细胞的信号流模式
p1 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
p2 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
p1 + p2


# 合并疾病和正常cellchat对象 ---------------------------------------------------
load("data/cellchat_Control.rdata")
object.list = list(Control = control.cellchat, Case = case.cellchat)
cellchat = mergeCellChat(object.list, add.names = names(object.list))

gg1 = compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 = compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
plot = gg1+gg2
plot
pdf("graph/case_control.pdf", width = 10, height = 6)
p1 = netVisual_bubble(cellchat,
                      sources.use = c("OCs"),
                      targets.use = c("Myeloid cells","Osteoblastic OS cells","NK/T cells","CAFs"),
                      comparison = c(1,2),
                      signaling = "SPP1",
                      remove.isolate = FALSE,
                      font.size = 14)
p1
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2
