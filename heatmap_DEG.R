library(Seurat)
library(ComplexHeatmap)
library(dplyr)
library(circlize)

dataRNA <- readRDS("scRNA_seq.rds")
dataRNA <- NormalizeData(dataRNA)
dataRNA <- FindVariableFeatures(dataRNA, nfeatures = 2000)
dataRNA <- ScaleData(dataRNA)
dataRNA <- RunPCA(dataRNA, npcs = 100)

scaled_matrix <- GetAssayData(dataRNA, layer = "scale.data")

rna_groups <- dataRNA@meta.data$class

df_up <- read.csv("UP.csv", header = TRUE)
df_down <- read.csv("DOWN.csv", header = TRUE)
outputFile <- "Heatmap.pdf"
figureTitle <- "ND vs T2D"

if (!"avg_log2FC" %in% colnames(df_up) | !"avg_log2FC" %in% colnames(df_down)) {
  stop("Error: avg_log2FC column not found in input files.")
}

df_up <- df_up %>% arrange(desc(avg_log2FC))  
df_down <- df_down %>% arrange(avg_log2FC)    

list1 <- df_up$X
list2 <- df_down$X

gene_cluster_info <- data.frame(
  Gene = c(list1, list2),
  Cluster = c(rep("UP", length(list1)), rep("DOWN", length(list2))),
  avg_log2FC = c(df_up$avg_log2FC, df_down$avg_log2FC)  
)

filtered_matrix <- scaled_matrix[rownames(scaled_matrix) %in% gene_cluster_info$Gene, ]

gene_cluster_info <- gene_cluster_info %>% 
  filter(Gene %in% rownames(filtered_matrix)) %>%
  arrange(desc(Cluster), desc(avg_log2FC))  # 'UP' first, sorted by avg_log2FC descending

filtered_matrix <- filtered_matrix[gene_cluster_info$Gene, ]

group_means <- aggregate(t(filtered_matrix), by = list(Class = rna_groups), FUN = mean)
rownames(group_means) <- group_means$Class
group_means <- group_means[, -1]

group_means["ND", ] <- colMeans(group_means[c("LN-ND", "OW-ND", "OB-ND"), ], na.rm = TRUE)
selected_means <- group_means[c("ND", "T2D"), ]

zscore_matrix <- t(scale(t(selected_means)))
zscore_matrix_transposed <- t(zscore_matrix)

col_fun <- colorRamp2(c(min(zscore_matrix_transposed), 0, max(zscore_matrix_transposed)), 
                      c("#0A5EB0", "#FDF7F4", "#F93827"))  

row_annotation <- rowAnnotation(
  Regulation = factor(gene_cluster_info$Cluster, levels = c("UP", "DOWN")),
  col = list(Regulation = c("UP" = "#A8CD89", "DOWN" = "#FADA7A")),  
  annotation_legend_param = list(title = "Regulation")
)

pdf(file = outputFile, width = 10, height = 12)

Heatmap(
  zscore_matrix_transposed,
  name = "Z-score",
  col = col_fun,  
  cluster_rows = FALSE,  
  cluster_columns = FALSE,  
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  row_title = "Genes",
  column_title = figureTitle,
  left_annotation = row_annotation 
)

dev.off()

