library(CCA)
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

commonID <- intersect(data1$GeneID, data2$GeneID)
data1 <- data1[data1$GeneID %in% commonID, ]
data2 <- data2[data2$GeneID %in% commonID, ]
data2 <- data2[match(data1$GeneID, data2$GeneID), ]

data1 <- data1[, !names(data1) %in% c("GeneID")]
data2 <- data2[, !names(data2) %in% c("GeneID")]

data1 <- as.data.frame(lapply(data1, as.numeric))
data2 <- as.data.frame(lapply(data2, as.numeric))

ccaResult <- cc(data1, data2)

xScores <- ccaResult$scores$xscores
yScores <- ccaResult$scores$yscores

xScoresDF <- data.frame(xScores)
yScoresDF <- data.frame(yScores)

mergedDF <- data.frame(
  CanonicalVariate1 = c(xScoresDF[, 1], yScoresDF[, 1]),
  CanonicalVariate2 = c(xScoresDF[, 2], yScoresDF[, 2]),
  Source = factor(c(rep("Gandal", nrow(xScoresDF)), rep("Sundari", nrow(yScoresDF))))
)

plot <- ggplot(mergedDF, aes(x = CanonicalVariate1, y = CanonicalVariate2, color = Source)) +
  geom_point(size = 3) +
  labs(title = "CCA", x = "Canonical Variate 1", y = "Canonical Variate 2") +
  scale_color_manual(values = c("Gandal" = "#1C1678", "Sundari" = "#C40C0C"))

pdf("CCACanonicalVariate.pdf")
print(plot)
dev.off()

canonicalCorrelations <- ccaResult$cor

pdf("CanonicalCorrelations.pdf")
barplot(canonicalCorrelations, 
        main = "Canonical Correlations", 
        xlab = "Canonical Dimensions", 
        ylab = "Correlation", 
        col = "#4793AF",
        ylim = c(0, 0.5))
dev.off()

