library(UpSetR)

inputFile1 <- "cluster_1_DOWN.csv"
inputFile2 <- "cluster_2_DOWN.csv"
inputFile3 <- "cluster_3_DOWN.csv"
inputFile4 <- "cluster_1_DOWN.csv"
inputFile5 <- "cluster_2_DOWN.csv"
inputFile6 <- "cluster_3_DOWN.csv"

inputFile7 <- "cluster_1_UP.csv"
inputFile8 <- "cluster_2_UP.csv"
inputFile9 <- "cluster_3_UP.csv"
inputFile10 <- "cluster_1_UP.csv"
inputFile11 <- "cluster_2_UP.csv"
inputFile12 <- "cluster_3_UP.csv"

outputPlot <- "upset_plot.pdf"

dataDF1 <- read.csv(inputFile1, header = TRUE, stringsAsFactors = FALSE)
dataDF2 <- read.csv(inputFile2, header = TRUE, stringsAsFactors = FALSE)
dataDF3 <- read.csv(inputFile3, header = TRUE, stringsAsFactors = FALSE)
dataDF4 <- read.csv(inputFile4, header = TRUE, stringsAsFactors = FALSE)
dataDF5 <- read.csv(inputFile5, header = TRUE, stringsAsFactors = FALSE)
dataDF6 <- read.csv(inputFile6, header = TRUE, stringsAsFactors = FALSE)
dataDF7 <- read.csv(inputFile7, header = TRUE, stringsAsFactors = FALSE)
dataDF8 <- read.csv(inputFile8, header = TRUE, stringsAsFactors = FALSE)
dataDF9 <- read.csv(inputFile9, header = TRUE, stringsAsFactors = FALSE)
dataDF10 <- read.csv(inputFile10, header = TRUE, stringsAsFactors = FALSE)
dataDF11 <- read.csv(inputFile11, header = TRUE, stringsAsFactors = FALSE)
dataDF12 <- read.csv(inputFile12, header = TRUE, stringsAsFactors = FALSE)

geneSets <- list(
  "Cluster 1 LN-ND vs OW-ND" = dataDF1$X,
  "Cluster 2 LN-ND vs OW-ND" = dataDF2$X,
  "Cluster 3 LN-ND vs OW-ND" = dataDF3$X,
  "Cluster 1 LN-ND vs OB-ND" = dataDF4$X,
  "Cluster 2 LN-ND vs OB-ND" = dataDF5$X,
  "Cluster 3 LN-ND vs OB-ND" = dataDF6$X,
  "Cluster 1 LN-ND vs T2D" = dataDF7$X,
  "Cluster 2 LN-ND vs T2D" = dataDF8$X,
  "Cluster 3 LN-ND vs T2D" = dataDF9$X,
  "Cluster 1 ND vs T2D" = dataDF10$X,
  "Cluster 2 ND vs T2D" = dataDF11$X,
  "Cluster 3 ND vs T2D" = dataDF12$X
)

binaryMatrix <- fromList(geneSets)

pdf(outputPlot, onefile = TRUE)
upset(binaryMatrix, 
      nsets = length(geneSets), 
      order.by = "freq", 
      main.bar.color = "#155E95", 
      sets.bar.color = "#bb0c00",
      keep.order = TRUE)
dev.off()

