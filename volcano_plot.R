library(ggplot2)
library(ggrepel)

inputFile <- "inputfile.csv"
outputFile <- "Volcano_Plot.pdf"
figureTitle <- "Volcano Plot"

volcanoData <- read.csv(inputFile, header = TRUE, stringsAsFactors = FALSE)

volcanoData$p_val_adj[volcanoData$p_val_adj == 0] <- 1e-300

volcanoData$expression <- "Not significant"
volcanoData$expression[volcanoData$avg_log2FC >= log2(1.5) & volcanoData$p_val_adj <= 0.1] <- "Upregulated"
volcanoData$expression[volcanoData$avg_log2FC <= -log2(1.5) & volcanoData$p_val_adj <= 0.1] <- "Downregulated"
volcanoData$expression[abs(volcanoData$avg_log2FC) < log2(1.5) & volcanoData$p_val_adj <= 0.1] <- "Significant but out of FC"

upregulated_genes <- subset(volcanoData, expression == "Upregulated")
downregulated_genes <- subset(volcanoData, expression == "Downregulated")

upregulatedCount <- nrow(upregulated_genes)
downregulatedCount <- nrow(downregulated_genes)
outOfFCCount <- nrow(subset(volcanoData, expression == "Significant but out of FC"))

legend_labels <- c(
  paste0("Upregulated (", upregulatedCount, ")"),
  paste0("Significant but out of FC (", outOfFCCount, ")"),
  paste0("Downregulated (", downregulatedCount, ")"),
  "Not significant"
)

xMin <- min(volcanoData$avg_log2FC)
xMax <- max(volcanoData$avg_log2FC)
yMin <- min(-log10(volcanoData$p_val_adj))
yMax <- max(-log10(volcanoData$p_val_adj)) + 1

pdf(file = outputFile, width = 10, height = 8)
ggplot(data = volcanoData, aes(x = avg_log2FC, y = -log10(p_val_adj), color = expression)) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), col = "#686D76", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.1), col = "#686D76", linetype = 'dashed') +
  geom_point(size = 1.5) +  # Use one geom_point for all data
  geom_text_repel(data = upregulated_genes, aes(label = X), size = 3, color = "#bb0c00", max.overlaps = 10) +
  geom_text_repel(data = downregulated_genes, aes(label = X), size = 3, color = "#155E95", max.overlaps = 10) +
  scale_color_manual(
    values = c(
      "Upregulated" = "#bb0c00",
      "Significant but out of FC" = "#050506",
      "Downregulated" = "#155E95",
      "Not significant" = "#B7B7B7"
    ),
    breaks = c("Upregulated", "Significant but out of FC", "Downregulated", "Not significant"),
    labels = legend_labels
  ) +
  coord_cartesian(ylim = c(yMin, yMax), xlim = c(xMin, xMax)) +
  labs(
    color = 'Expression',
    x = expression("log"[2]*"(FC)"),
    y = expression("-log"[10]*"(adj pValue)")
  ) +
  ggtitle(figureTitle)
dev.off()

