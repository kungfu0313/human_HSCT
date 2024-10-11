library(limma)
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

##-------- get differentially expressed genes of each cluster 

setwd("/PATH/TO/RDS/")
scRNA <- readRDS("./pre-processed.rds")


diff.wilcox = FindAllMarkers(scRNA)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "top10_diff_genes_wilcox.csv", row.names = F)

##-------- Celltyping by singleR
#hpca.se_fine
load("DB/hpca_se.RData")
 testdata <- GetAssayData(scRNA, slot="data")
 refdata <- hpca.se
 clusters <- scRNA@meta.data$seurat_clusters
 cellpred <- SingleR(test=testdata, ref=refdata, labels=refdata$label.fine, method="cluster", clusters=clusters, assay.type.test="logcounts", assay.type.ref="logcounts")
 celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
 write.csv(celltype,"celltype_singleR_hpca_ca.csv",row.names = F)
 scRNA@meta.data$celltypehpca = "NA"
 for(i in 1:nrow(celltype)){ scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltypehpca'] <- celltype$celltype[i]}
#hpca.se_main
load("DB/hpca_se.RData")
 testdata <- GetAssayData(scRNA, slot="data")
 refdata <- hpca.se
 clusters <- scRNA@meta.data$seurat_clusters
 cellpred <- SingleR(test=testdata, ref=refdata, labels=refdata$label.main, method="cluster", clusters=clusters, assay.type.test="logcounts", assay.type.ref="logcounts")
 celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
 write.csv(celltype,"celltype_singleR_hpca_main_ca.csv",row.names = F)
 scRNA@meta.data$celltypehpca_main = "NA"
 for(i in 1:nrow(celltype)){ scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltypehpca_main'] <- celltype$celltype[i]}
#blue.se_fine
load("DB/blue.se.RData")
 testdata <- GetAssayData(scRNA, slot="data")
 refdata <- Blue.se
 clusters <- scRNA@meta.data$seurat_clusters
 cellpred <- SingleR(test=testdata, ref=refdata, labels=refdata$label.fine, method="cluster", clusters=clusters, assay.type.test="logcounts", assay.type.ref="logcounts")
 celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
 write.csv(celltype,"celltype_singleR_blue_ca.csv",row.names = F)
 scRNA@meta.data$celltypeblue = "NA"
 for(i in 1:nrow(celltype)){ scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltypeblue'] <- celltype$celltype[i]}
#blue.se_main
load("DB/blue.se.RData")
 testdata <- GetAssayData(scRNA, slot="data")
 refdata <- Blue.se
 clusters <- scRNA@meta.data$seurat_clusters
 cellpred <- SingleR(test=testdata, ref=refdata, labels=refdata$label.main, method="cluster", clusters=clusters, assay.type.test="logcounts", assay.type.ref="logcounts")
 celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
 write.csv(celltype,"celltype_singleR_blue_main_ca.csv",row.names = F)
 scRNA@meta.data$celltypeblue_main = "NA"
 for(i in 1:nrow(celltype)){ scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltypeblue_main'] <- celltype$celltype[i]}
#DB_DBice.RData
load("DB/DB_DBice.RData")
 testdata <- GetAssayData(scRNA, slot="data")
 refdata <- DBice
 clusters <- scRNA@meta.data$seurat_clusters
 cellpred <- SingleR(test=testdata, ref=refdata, labels=refdata$label.main, method="cluster", clusters=clusters, assay.type.test="logcounts", assay.type.ref="logcounts")
 celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
 write.csv(celltype,"celltype_singleR_DBice.csv",row.names = F)
 scRNA@meta.data$celltypeDBice = "NA"
 for(i in 1:nrow(celltype)){ scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltypeDBice'] <- celltype$celltype[i]}
#DB_Nover.RData
load("DB/DB_Nover.RData")
 testdata <- GetAssayData(scRNA, slot="data")
 refdata <- Nover
 clusters <- scRNA@meta.data$seurat_clusters
 cellpred <- SingleR(test=testdata, ref=refdata, labels=refdata$label.main, method="cluster", clusters=clusters, assay.type.test="logcounts", assay.type.ref="logcounts")
 celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
 write.csv(celltype,"celltype_singleR_Nover.csv",row.names = F)
 scRNA@meta.data$celltypeNover = "NA"
 for(i in 1:nrow(celltype)){ scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltypeNover'] <- celltype$celltype[i]}
#DB_Monaco.RData
load("DB/DB_Monaco.RData")
 testdata <- GetAssayData(scRNA, slot="data")
 refdata <- Monaco
 clusters <- scRNA@meta.data$seurat_clusters
 cellpred <- SingleR(test=testdata, ref=refdata, labels=refdata$label.main, method="cluster", clusters=clusters, assay.type.test="logcounts", assay.type.ref="logcounts")
 celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
 write.csv(celltype,"celltype_singleR_Monaco.csv",row.names = F)
 scRNA@meta.data$celltypeMonaco = "NA"
 for(i in 1:nrow(celltype)){ scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltypeMonaco'] <- celltype$celltype[i]}


##-------- update metadata by adding cell type info of each cluster by combining both celltype markers and singleR predictions

scRNA@meta.data$celltype <- 'NA'
for(i in 1:length(scRNA@meta.data$celltype))
{if(scRNA@meta.data$seurat_clusters[i] == 0) scRNA@meta.data$celltype[i] <- 'NK'
else if (scRNA@meta.data$seurat_clusters[i] == 1) scRNA@meta.data$celltype[i] <- 'T cell'
else if (scRNA@meta.data$seurat_clusters[i] == 2) scRNA@meta.data$celltype[i] <- 'Pre B cell'
else if (scRNA@meta.data$seurat_clusters[i] == 3) scRNA@meta.data$celltype[i] <- 'Macrophage'
else if (scRNA@meta.data$seurat_clusters[i] == 4) scRNA@meta.data$celltype[i] <-  'Monocyte'
...
}

##-------- save celltyping results  
saveRDS(scRNA, "./celltyped.rds")
