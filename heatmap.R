library(pheatmap)
library(ggplot2)


##-------- Generating heatmaps( Figure 5,  Figure S16, Figure S17, and etc.)
 

dataset <- read.table(file="input",sep="\t", header = T, row.names=1,check.names=FALSE)  #all input files in this study were pre-processed, they were either ratio or pre-scaled.
exp_ds = dataset[, c(1:ncol)]
p <- pheatmap(exp_ds, show_rownames=T, show_colnames=T, diBMay_numbers=F,color=colorRampPalette(c("navy", "white","firebrick3"))(50), cluster_col = F,cluster_rows=F, main="heatmap", cellheight=27, cellwidth=27, fontsize=20)
ggsave("heatmap.pdf", p, width=25,height=20)
ggsave("heatmap.bmp", p, width=25,height=20)
