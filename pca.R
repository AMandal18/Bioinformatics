library(ggplot2)
library(reshape2)

data1 <- read.table("data1.txt", header = TRUE)
data1$GeneID <- sub("\\..*", "", data1$GeneID)
data1 <- data1[ , !(colnames(data1) == "GeneName")]
data1 <- as.data.frame(data1)

data2 <- read.table("data2.txt", header = TRUE)
data2$GeneID <- sub("\\..*", "", data2$GeneID)
data2 <- data2[ , !(colnames(data2) == "GeneName")]
data2 <- as.data.frame(data2)

mergedData <- merge(data1, data2, by = "GeneID")
mergedDataNew <- mergedData[, !names(mergedData) %in% c("GeneID")]
numericData <- as.data.frame(lapply(mergedDataNew, as.numeric))

numericDataT <- t(numericData)
pcaResult <- prcomp(numericDataT, scale. = TRUE)

pcScores <- pcaResult$x

mergedDF <- data.frame(
  PC1 = pcScores[, 1],
  PC2 = pcScores[, 2],
  Source = factor(c(rep("Gandal", 808), rep("Sundari", 856 - 809 + 1)))
)

plot <- ggplot(mergedDF, aes(x = PC1, y = PC2, color = Source)) +
  geom_point(size = 3) +
  labs(title = "PCA", x = "Principal Component 1", y = "Principal Component 2") +
  scale_color_manual(values = c("Gandal" = "#1C1678", "Sundari" = "#C40C0C"))

pdf("/data/thakurelaLab/users/Ankita/Sundari/Analysis/PCA.pdf")
print(plot)
dev.off()

