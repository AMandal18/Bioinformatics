library(ggplot2)
library(dplyr)
library(reshape2)

inputFile1 <- "Cluster_1_UP.csv"
inputFile2 <- "Cluster_1_DOWN.csv"
inputFile3 <- "Cluster_2_UP.csv"
inputFile4 <- "Cluster_2_DOWN.csv"
inputFile5 <- "Cluster_3_UP.csv"
inputFile6 <- "Cluster_3_DOWN.csv"

outputFile <- "Bubble_Plot.pdf"

dataDF1 <- read.csv(inputFile1, header = TRUE, stringsAsFactors = FALSE)
dataDF2 <- read.csv(inputFile2, header = TRUE, stringsAsFactors = FALSE)
dataDF3 <- read.csv(inputFile3, header = TRUE, stringsAsFactors = FALSE)
dataDF4 <- read.csv(inputFile4, header = TRUE, stringsAsFactors = FALSE)
dataDF5 <- read.csv(inputFile5, header = TRUE, stringsAsFactors = FALSE)
dataDF6 <- read.csv(inputFile6, header = TRUE, stringsAsFactors = FALSE)

allData <- list(dataDF1, dataDF2, dataDF3, dataDF4, dataDF5, dataDF6)

labels <- c("Cluster 1 UP", "Cluster 1 DOWN", "Cluster 2 UP", "Cluster 2 DOWN", "Cluster 3 UP", "Cluster 3 DOWN")

allData <- lapply(seq_along(allData), function(i) {
  data <- allData[[i]]
  data$Dataset <- labels[i]
  return(data)
})

combinedData <- do.call(rbind, allData)

combinedData$Dataset <- factor(combinedData$Dataset, levels = labels)

bubblePlot <- ggplot(combinedData, aes(x = Dataset, y = X)) +
  geom_point(aes(size = avg_log2FC, color = p_val_adj)) +
  scale_color_gradient2(low = "#A294F9", mid = "#F0BB78", high = "#FF8383", midpoint = 0) +
  scale_size_continuous(range = c(3, 10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "LN-ND vs T2D",
    x = "Condition",
    y = "Gene",
    size = "Log Fold Change",
    color = "p-value"
  )

ggsave(outputFile, bubblePlot, width = 12, height = 22)

