library(Seurat)
library(sankeywheel)
library(patchwork)
library(limma)
library(SingleR)
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)
library(stringr)
library(DT)
library(sankeywheel)
library(hdf5r)


##-------- Generating figures( Figure 2B,  Figure 2C, Figure 4, and etc.)

setwd("/PATH/TO/RDS/")
scRNA <- readRDS("./celltyped.rds")

##-------- figure by celltype (TSNE)
p1<-DimPlot(scRNA, reduction = "tsne", label = TRUE, group.by ="celltype", pt.size = 1.5, label.size =15,repel=T)
ggsave("RGD007_tsne_bycelltype.pdf", p1, width=23, height=17)	
p1<-DimPlot(scRNA, reduction = "tsne", label = FALSE, group.by ="celltype", pt.size = 1.5, label.size =15,  repel=T)
ggsave("tsne_bycelltype_nolabel.pdf", p1, width=23,height=17)


##-------- Dotplot for markers 
new.cluster.ids <- factor(scRNA@meta.data$celltype, levels=unique(scRNA@meta.data$celltype))
names(new.cluster.ids) <- rownames(scRNA@meta.data)
scRNA@active.ident <- new.cluster.ids

#Please  update level list based on real celltype list
my_level <- c('CD4+ T', 'CD8+ T', 'gd T', 'Treg', 'Macrophage' ,'Monocyte' ,'Basophils' ,'Granulocytes' ,'Mast cell' ,'B cell','NK') 
Idents(scRNA) <- factor(Idents(scRNA), levels=my_level)

#Please update celltype marker list based on real celltype list
select_genes_all=c('FCGR3A','IFNG','FCMR','CD79A','CLC','OSM','LY86','C1orf162','GATA2','TPSAB1','SYAP1','TEC','CD14','S100A12','CTLA4','FOXP3','CCR9','TRDC', 'RPS8','TNF','TNFRSF4','ICOS','CD33','CD19','CD8B','CD8A','CD4','CD3G','PTPRC')

p<-DotPlot(scRNA,features=select_genes_all, assay="RNA", col.min=-1) + theme(axis.text.x=element_text(angle=45,vjust=0.95, hjust=0.95,size=18),axis.text.y=element_text(size=18) ) + RotatedAxis() + coord_flip()+scale_color_gradient2(high="darkred",mid = "pink",low ="lightblue", midpoint = 0)

ggsave("dotplot_markers.pdf", p, width=8 ,height=15)

