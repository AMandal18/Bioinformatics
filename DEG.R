library(Seurat)
library(dplyr)
library(ggplot2)

dataRNA <- readRDS("scRNA_seq.rds")

dataRNA <- NormalizeData(dataRNA)
dataRNA <- FindVariableFeatures(dataRNA, nfeatures = 2000)
dataRNA <- ScaleData(dataRNA)
dataRNA <- RunPCA(dataRNA, npcs = 100)

scaled_matrix <- GetAssayData(dataRNA, layer = "scale.data")

rna_groups <- dataRNA@meta.data$class

list1 <- read.csv("Cluster1_Genes.csv", header = FALSE, stringsAsFactors = FALSE)[,1]
list2 <- read.csv("Cluster2_Genes.csv", header = FALSE, stringsAsFactors = FALSE)[,1]
list3 <- read.csv("Cluster3_Genes.csv", header = FALSE, stringsAsFactors = FALSE)[,1]

Idents(dataRNA) <- "class"

dge_cluster1_LNOW <- FindMarkers(
  object = dataRNA,
  ident.1 = "OW-ND",
  ident.2 = "LN-ND", 
  features = list1, 
  test.use = "wilcox",
  min.pct = 0.1,  
  logfc.threshold = 0 
)
dge_cluster1_LNOW$p_val_adj <- p.adjust(dge_cluster1_LNOW$p_val, method = "BH")
dge_cluster1_LNOW$expression <- NA
dge_cluster1_LNOW$expression[dge_cluster1_LNOW$avg_log2FC >= log2(1.5) & dge_cluster1_LNOW$p_val_adj <= 0.1] <- "Upregulated"
dge_cluster1_LNOW$expression[dge_cluster1_LNOW$avg_log2FC <= -log2(1.5) & dge_cluster1_LNOW$p_val_adj <= 0.1] <- "Downregulated"
write.csv(dge_cluster1_LNOW, file = "cluster_1.csv", row.names = TRUE)
upregulated_genes <- subset(dge_cluster1_LNOW, expression == "Upregulated")
write.csv(upregulated_genes, file = "cluster_1_Upregulated.csv", row.names = TRUE)
downregulated_genes <- subset(dge_cluster1_LNOW, expression == "Downregulated")
write.csv(downregulated_genes, file = "cluster_1_Downregulated.csv", row.names = TRUE)

