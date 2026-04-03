# ===== 镜像设置（国内服务器使用）=====
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/BiocMirror/")

# ===== 并行计算设置 =====
library(future)
# 设置并行策略为 multisession，核数为 4
plan("multisession", workers = 4)
# 允许传输较大的对象（单细胞数据通常很大，建议设置为 50GB 以上）
options(future.globals.maxSize = 50 * 1024^3)

# 安装所需的包（如需要）
# install.packages(c("dplyr", "data.table", "ggplot2", "Seurat", "tidyverse", "clustree", "patchwork", "harmony", "GO.db", "clusterProfiler", "org.Hs.eg.db", "GOplot", "ggpubr", "ggstar"))

library(dplyr)
library(data.table)
library(ggplot2)
library(Seurat)
library(tidyverse)
library(clustree)
# library(decontX)
library(patchwork)
library(harmony)

rm(list = ls())
gc()
setwd("~/autodl-tmp/Single _cell")

load("seurat_object.rdata")
# 差异分析 ---------------------------------------------------------------------
Idents(seurat_object)
Idents(seurat_object) = seurat_object$cell.types  #用Idents定位到细胞类型
# 组间差异分析
seurat_object$celltype.group = paste(seurat_object$cell.types, seurat_object$condition, sep = "_")  #细胞类型和分组信息合并
# 设置Idents为celltype.group，以便差异分析使用
Idents(seurat_object) = seurat_object$celltype.group
head(seurat_object$celltype.group)
table(seurat_object$cell.types, seurat_object$condition)
table(seurat_object$celltype.group)
# save(seurat_object, file = "scRNA.rdata")

# 先对每个细胞类型分别设置Idents，然后做差异分析
# CAFs 差异分析
Idents(seurat_object) = seurat_object$celltype.group
marker1 = FindMarkers(seurat_object,
                      ident.1 = "CAFs_EarlyAF",      #指明做哪个细胞类型的疾病和正常组分析
                      ident.2 = "CAFs_PermanentAF",
                      logfc.threshold = 0.1,  #差异倍数最小值
                      min.pct = 0.1  #设置在两个群体中任一群体的最少比例细胞中检测到的基因的阈值
                      )
write.csv(marker1, "data/CAFs_diff.csv")

# 循环跑差异分析 - 对每个细胞类型
cellfordeg = unique(seurat_object$cell.types)
for(i in 1:length(cellfordeg)){
  print(i)
  cellDEG = FindMarkers(seurat_object,
                        ident.1 = paste0(cellfordeg[i],"_EarlyAF"),
                        ident.2 = paste0(cellfordeg[i],"_PermanentAF"),
                        logfc.threshold = 0.1,
                        min.pct = 0.1)
  # 处理文件名中的特殊字符
  cellname = gsub("[/ ]", "_", cellfordeg[i])
  write.csv(cellDEG, paste0("data/", cellname, "_diff.csv"))
}

# 对指定的俩个聚类细胞进行差异分析
# 先设置正确的Idents
Idents(seurat_object) = seurat_object$cell.types
marker2 = FindMarkers(seurat_object,
                      ident.1 = "Osteoclasts OCs",
                      ident.2 = "CAFs",
                      logfc.threshold = 0.1,
                      min.pct = 0.1)
write.csv(marker2, "data/OsteoclastsOCs_CAFs_diff.csv")

# 绘制火山图
library(ggplot2)
pval = 0.05
logFC = 1
diff = marker1
diff$change = ifelse((diff$p_val_adj < pval & abs(diff$avg_log2FC) > logFC),
                     ifelse(diff$avg_log2FC > logFC, "Up", "Down"), "None")
table(diff$change)
write.csv(diff, "data/CAFs_DEGs.csv")

p1 = ggplot(data = diff, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = change)) + # P值取负对数,相当于越显著,负对数越大
  geom_point(alpha = 1) +  # 点的透明度
  scale_color_manual(values = c("#3288BD", "#808080","#D53E4F")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 12)) + 
  ylab(expression(-log[10]("P.Value"))) +    # expression的作用就是让log10的10下标
  xlab(expression(log[2]("Fold Change"))) +
  geom_vline(xintercept = c(-1, 1),          # 加垂直线
             lty = 2,                        # 线的类型
             col = "black",                  # 颜色
             lwd = 0.7) +                    # 宽度
  geom_hline(yintercept = -log10(0.05),
             lty = 2,
             col = "black",
             lwd = 0.7)
p1
ggsave("graph/Volcano_plot.pdf", width = 8, height = 6, plot = p1)
ggsave("graph/Volcano_plot.tiff", width = 8, height = 6, plot = p1, dpi = 300)


## 富集分析
library(GO.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)
library(GOplot)
library(org.Mm.eg.db)
library(ggpubr)
library(ggstar)
# 富集分析
DEGs = read.csv("./data/CAFs_DEGs.csv", row.names = 1, check.names = F)
DEGs_up = DEGs[DEGs$change %in% "Up",]
DEGs_down = DEGs[DEGs$change %in% "Down",]

gene_id = bitr(rownames(DEGs_up),         #输入文件
               fromType = "SYMBOL",       #转换前的基因ID类型
               toType =  c("ENTREZID"),   #转换后的基因ID类型
               OrgDb = org.Hs.eg.db)

#GO富集分析
GO_all = enrichGO(gene_id$ENTREZID,        #输入基因的"ENTREZID"
                  OrgDb = org.Hs.eg.db,    #注释信息
                  keyType = "ENTREZID",
                  ont = "all",             #可选条目BP,CC,MF
                  pAdjustMethod = "BH",    #p值的校正方式
                  pvalueCutoff = 1,     #pvalue的阈值
                  qvalueCutoff = 1,     #qvalue的阈值
                  readable = TRUE)         #是否将ENTREZID转换为symbol
GO_all_result = data.frame(GO_all)
GO_all_result = GO_all_result[GO_all_result$pvalue < 0.05,]
write.csv(GO_all_result, "data/GO_result_up.csv")

# GO可视化
GO = read.csv("data/GO_result_up.csv", row.names = 1, check.names = F)
GO=GO[1:16,]
GO$Description = factor(GO$Description, levels = GO$Description)
p1 = ggplot(GO, aes(x = Description, y = Count, color = pvalue))+
  geom_segment(aes(x = Description, xend = Description, y = 0, yend = Count), size = 1.8, color = "gray30")+
  geom_star(aes(fill = pvalue), size = 8)+
  scale_color_gradient(low = "#f1aea7", high = "#abc8e5")+
  scale_fill_gradient(low = "#f1aea7", high = "#abc8e5")+
  labs(title = "The Most Enriched GO Terms",     #设置标题、x轴和Y轴名称
       x = "", 
       y = "Gene number")+
  theme(axis.text.y = element_text(size = 17, colour = c(rep("#c25759", 6),rep("#599cb4", 6))),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 18, colour = "grey20"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length = unit(1.5, "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.background = element_blank(),
        panel.grid.major = element_line(linewidth = 0.8, color = "#E8E8E8"),
        plot.background = element_rect(color = "white"),
        plot.title = element_text(size = 16))+
  coord_flip()
p1
ggsave("graph/GO_result.pdf", width = 8.6, height = 5, plot = p1)
ggsave("graph/GO_result.tiff", width = 8.6, height = 5, dpi = 300, plot = p1)

# KEGG富集分析
KEGG = enrichKEGG(gene_id$ENTREZID,
                  keyType = "kegg",
                  organism = "hsa",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1,
                  pAdjustMethod = "BH")
KEGG_result = setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") %>% data.frame()
KEGG_result = KEGG_result[KEGG_result$pvalue < 0.05,]
write.csv(KEGG_result, "data/KEGG_result_up.csv")

# KEGG可视化
KEGG = read.csv("GOKEGG/KEGG.csv", row.names = 1, check.names = F)
KEGG$Description = factor(KEGG$Description, levels = KEGG$Description)
p2 = ggplot(KEGG, aes(x = Description, y = Count, color = pvalue))+
  geom_segment(aes(x = Description, xend = Description, y = 0, yend = Count), size = 1.8, color = "gray30")+
  geom_star(aes(fill = pvalue), size = 8)+
  scale_color_gradient(low = "#f1aea7", high = "#abc8e5")+
  scale_fill_gradient(low = "#f1aea7", high = "#abc8e5")+
  labs(title = "The Most Enriched KEGG Terms",     #设置标题、x轴和Y轴名称
       x = "", 
       y = "Gene number")+
  theme(axis.text.y = element_text(size = 17, colour = c(rep("#c25759", 6),rep("#599cb4", 6))),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 18, colour = "grey20"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length = unit(1.5, "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.background = element_blank(),
        panel.grid.major = element_line(linewidth = 0.8, color = "#E8E8E8"),
        plot.background = element_rect(color = "white"),
        plot.title = element_text(size = 16))+
  coord_flip()
p2
ggsave("graph/KEGG_result.pdf", width = 8.6, height = 5, plot = p2)
ggsave("graph/KEGG_result.tiff", width = 8.6, height = 5, dpi = 300, plot = p2)



