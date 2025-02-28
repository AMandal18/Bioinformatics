library(DOSE)           
library(clusterProfiler) 
library(org.Hs.eg.db)    
library(enrichplot)      
library(ggplot2)         
library(data.table)      

inputFileDF <- "inputfile.csv"
outputFileDot <- "DO_dot.pdf"
outputFileBar <- "DO_bar.pdf"
outputFileTable <- "DO_results.csv"

dataDF <- fread(inputFileDF)

gene_list <- dataDF$X

gene_entrez <- mapIds(org.Hs.eg.db, 
                      keys = gene_list, 
                      column = "ENTREZID", 
                      keytype = "SYMBOL", 
                      multiVals = "first")

gene_entrez <- na.omit(gene_entrez)

do_result <- enrichDO(gene = gene_entrez, 
                      pvalueCutoff = 0.05)

print(head(do_result))

write.csv(do_result@result, file = outputFileTable, row.names = FALSE)

pdf(outputFileDot, width = 8, height = 6)
dotplot(do_result, showCategory=10)
dev.off()

pdf(outputFileBar, width = 8, height = 6)
barplot(do_result, showCategory=10)
dev.off()

