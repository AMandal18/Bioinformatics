library(STRINGdb)
library(org.Hs.eg.db)  
library(ggplot2)      
library(clusterProfiler)

inputFileDF <- "inputfile.csv"
outputFilePPI <- "string_ppi.pdf"
outputFilePath <- "string_pathways.pdf"

dataDF <- fread(inputFileDF)

gene_list <- dataDF$X

string_db <- STRINGdb$new( version="12.0", species=9606, score_threshold=200, input_directory="")

gene_entrez <- mapIds(org.Hs.eg.db, 
                      keys = gene_list, 
                      column = "ENTREZID", 
                      keytype = "SYMBOL", 
                      multiVals = "first")

gene_entrez <- na.omit(gene_entrez)

genes_df <- data.frame(EntrezID = gene_entrez)

mapped_genes <- string_db$map(genes_df, "EntrezID", removeUnmappedRows=TRUE)
print(head(mapped_genes))

ppi_interactions <- string_db$get_interactions(mapped_genes$STRING_id)
print(head(ppi_interactions))

pdf(outputFilePPI, width = 10, height = 8)  
string_db$plot_network(mapped_genes$STRING_id)
dev.off()

enrichment_results <- string_db$get_enrichment(mapped_genes$STRING_id)
print(head(enrichment_results))

sig_pathways <- enrichment_results[enrichment_results$fdr < 0.05, ]

p <- ggplot(sig_pathways[1:20,], aes(x=reorder(description, -log10(fdr)), y=-log10(fdr))) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  labs(title="Top 20 Enriched STRING Pathways", x="Pathway", y="-log10(p-value)") +
  theme_minimal()
ggsave(outputFilePath, plot = p, width = 10, height = 6)

