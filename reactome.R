library(ReactomePA)      
library(clusterProfiler) 
library(org.Hs.eg.db)    
library(enrichplot)      
library(ggplot2)         
library(data.table)     

inputFileDF <- "inputfile.csv"
outputFileDot <- "Reactome_dot.pdf"
outputFileBar <- "Reactome_bar.pdf"
outputFileTable <- "Reactome_results.csv"

dataDF <- fread(inputFileDF)

gene_list <- dataDF$X

gene_entrez <- mapIds(org.Hs.eg.db, 
                      keys = gene_list, 
                      column = "ENTREZID", 
                      keytype = "SYMBOL", 
                      multiVals = "first")

gene_entrez <- na.omit(gene_entrez)

reactome_result <- enrichPathway(gene = gene_entrez, 
                                 organism = "human", 
                                 pvalueCutoff = 0.05)

print(head(reactome_result))

write.csv(reactome_result@result, file = outputFileTable, row.names = FALSE)

pdf(outputFileDot, width = 8, height = 6)
dotplot(reactome_result, showCategory=10)
dev.off()

pdf(outputFileBar, width = 8, height = 6)
barplot(reactome_result, showCategory=10)
dev.off()

