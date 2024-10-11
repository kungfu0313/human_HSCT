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

##-------- Read in cellranger filtered data 
setwd("/PATH/TO/CELLRANGER/OUTS/")
data <- Read10X(data.dir = "filtered_feature_bc_matrix")
scRNA = CreateSeuratObject(counts = data, min.cells = 0,  min.features = 200, project = "cellranger")

##-------- QC
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
#head(scRNA@meta.data)

#Vlnplot of percent.mt, nFeature_RNA and percent.HB before filtering
col.num <- length(levels(scRNA@active.ident))
violin <- VlnPlot(scRNA,
 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
 cols =rainbow(col.num),
 pt.size = 0.01,
 ncol = 4) + 
 theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
 ggsave("vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
 ggsave("vlnplot_before_qc.png", plot = violin, width = 12, height = 6)  
 plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
 plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
 plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
 pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
 ggsave("pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5) 
 ggsave("pearplot_before_qc.png", plot = pearplot, width = 12, height = 5)

#Vlnplot of percent.mt, nFeature_RNA and percent.HB after filtering
#minGene=200£¬maxGene=7500£¬pctMT=10¡£
minGene=200
maxGene=7500
pctMT=10
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(scRNA@active.ident))
violin <-VlnPlot(scRNA,
 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
 cols =rainbow(col.num), pt.size = 0.1,ncol = 4) + 
 theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
 ggsave("vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 
 ggsave("vlnplot_after_qc.png", plot = violin, width = 12, height = 6)


##--------checking Variable Features
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) 

##--------Scale data
scale.genes <-  rownames(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes) 


##--------Cell Cycle effect related
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA))

#Cell Cycle scoring 
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA))
scRNA <- CellCycleScoring(object=scRNA,  g2m.features=g2m_genes,  s.features=s_genes)

##Checking Cell Cycle effect
scRNAa <- RunPCA(scRNA, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
ggsave("cellcycle_pca.pdf", p, width = 8, height = 6)

##Removing Cell cycle effect by regressing out cell cycle scores if needed.
scRNA <- ScaleData(scRNA, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(scRNA))

#checking Cell Cycle effect again
scRNAa <- RunPCA(scRNA, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
ggsave("cellcycle_pca_afterregress.pdf", p, width = 8, height = 6)


##--------PCA
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(scRNA, ndims=50, reduction="pca")
plotc <- plot1+plot2
ggsave("pca.pdf", plot = plotc, width = 12, height = 5)

#pick proper PCA number based on pca.pdf, for example 40, and do cell clustering
pc.num=1:40
scRNA <- FindNeighbors(scRNA, dims = pc.num)
scRNA <- FindClusters(scRNA, resolution = 0.7)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- FindClusters(scRNA, resolution = 0.6)
table(scRNA@meta.data$seurat_clusters)

#-------- TSNE or UMAP
pc.num=1:40  #resolution: 0.6
#tSNE
scRNA = RunTSNE(scRNA, dims = pc.num)
#UMAP
scRNA <- RunUMAP(scRNA, dims = pc.num)

##-------- save pre-processed data as RDS 
saveRDS(scRNA, "./pre-processed.rds")
