library(Seurat)
library(ggplot2)
library(MASS) 

dataRNA <- readRDS("scRNA_seq.rds")

dataRNA <- NormalizeData(dataRNA)
dataRNA <- FindVariableFeatures(dataRNA, nfeatures = 2000)
dataRNA <- ScaleData(dataRNA)
dataRNA <- RunPCA(dataRNA, npcs = 100)

scaled_matrix <- GetAssayData(dataRNA, layer = "scale.data")

rna_groups <- dataRNA@meta.data$class

list1 <- read.csv("Cluster_1.csv", header = TRUE)$Gene
list2 <- read.csv("Cluster_2.csv", header = TRUE)$Gene
list3 <- read.csv("Cluster_3.csv", header = TRUE)$Gene

cluster_genes <- list(
  Cluster1 = list1,
  Cluster2 = list2,
  Cluster3 = list3
)

group_means <- lapply(cluster_genes, function(cluster_genes_list) {
  cluster_data <- subset(dataRNA, features = cluster_genes_list)
  
  if (length(rownames(cluster_data)) == 0) {
    warning("No genes found in cluster for the subset.")
    return(NULL)
  }
  
  cluster_expression <- GetAssayData(cluster_data, assay = "alpha_RNA", layer = "counts")
  
  aggregate_data <- aggregate(
    t(as.matrix(cluster_expression)),
    by = list(Group = dataRNA@meta.data$class),
    FUN = mean
  )
  
  rownames(aggregate_data) <- aggregate_data$Group
  aggregate_data$Group <- NULL
  return(aggregate_data)
})

group_means <- Filter(Negate(is.null), group_means)

all_genes <- Reduce(union, lapply(group_means, colnames))

aligned_group_means <- lapply(group_means, function(df) {
  missing_genes <- setdiff(all_genes, colnames(df))
  
  for (gene in missing_genes) {
    df[[gene]] <- NA
  }
  
  df <- df[, all_genes, drop = FALSE]
  
  return(df)
})

mean_expression <- do.call(rbind, aligned_group_means)

mean_expression[is.na(mean_expression)] <- 0

lda_input <- as.data.frame(mean_expression)
lda_input$Group <- rep(names(group_means), each = nrow(group_means[[1]]))

lda_result <- lda(Group ~ ., data = lda_input)

lda_scores <- predict(lda_result)$x

lda_df <- as.data.frame(lda_scores)

lda_df$Cluster <- rep(names(group_means), each = nrow(group_means[[1]]))
lda_df$Group <- sapply(rownames(lda_df), function(x) {
  strsplit(x, "\\.")[[1]][2]  
})

custom_colors <- c("LN-ND" = "#C62E2E", "OB-ND" = "#31511E", "OW-ND" = "#F09319", "T2D" = "#8B5DFF")

lda_plot <- ggplot(lda_df, aes(x = LD1, y = LD2, color = Group, shape = Cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "LDA of Mean Gene Expression",
    x = "LD1",
    y = "LD2"
  ) +
  scale_color_manual(values = custom_colors) +  
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

pdf("LDA.pdf", 
    width = 8, height = 6)
print(lda_plot)
dev.off()

