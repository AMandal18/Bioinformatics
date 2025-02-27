library(clusterProfiler)
library(org.Hs.eg.db)  
library(enrichplot)
library(data.table)

inputFileDF <- "inputfile.csv"
outputFileDot <- "KEGG_dot.pdf"
outputFileBar <- "KEGG_bar.pdf"
outputFileCnet <- "KEGG_cnet.pdf"
outputFileEmap <- "KEGG_emap.pdf"

dataDF <- fread(inputFileDF)

gene_list <- dataDF$X

gene_entrez <- mapIds(org.Hs.eg.db, 
                      keys = gene_list, 
                      column = "ENTREZID", 
                      keytype = "SYMBOL", 
                      multiVals = "first")

gene_entrez <- na.omit(gene_entrez)

gene_kegg <- bitr_kegg(gene_entrez, fromType = "ncbi-geneid", toType = "kegg", organism = "hsa")

kegg_result <- enrichKEGG(gene = gene_kegg$kegg, 
                          organism = "hsa", 
                          keyType = "kegg",
                          pvalueCutoff = 0.1, 
                          qvalueCutoff = 0.2)  

print(head(kegg_result))

pdf(outputFileDot, width = 8, height = 6)
dotplot(kegg_result, showCategory=10)
dev.off()

pdf(outputFileBar, width = 8, height = 6)
barplot(kegg_result, showCategory=10)
dev.off()

gene_symbol_map <- mapIds(org.Hs.eg.db, 
                          keys = gene_entrez, 
                          column = "SYMBOL", 
                          keytype = "ENTREZID", 
                          multiVals = "first")

kegg_result@result$geneID <- sapply(strsplit(kegg_result@result$geneID, "/"), function(x) {
  paste(na.omit(gene_symbol_map[x]), collapse = "/")
})

kegg_result <- pairwise_termsim(kegg_result)

pdf(outputFileCnet, width = 8, height = 6)  
cnetplot(kegg_result, showCategory = 10, node_label = "all")  
dev.off()

pdf(outputFileEmap, width = 8, height = 6)
emapplot(kegg_result)
dev.off()

