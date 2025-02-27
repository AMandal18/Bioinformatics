library(data.table)
library(topGO)

inputFileDF1 <- "GeneList.tsv"

inputFileDF2 <- "file.csv"
outputFileDF <- "TopGO.tsv"

geneNamesDF <- fread(inputFileDF1)
dataDF <- fread(inputFileDF2)

geneNames <- unique(sort(geneNamesDF[[1]]))
myInterestingGenes <- unique(sort(dataDF[[1]]))

geneList=factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList)=geneNames
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="SYMBOL")

GOdata = new('topGOdata',
         ontology = 'BP',
         allGenes = geneList,
         annot = annFUN.GO2genes,
         GO2genes=allGO2genes,
         geneSel = myInterestingGenes,
         nodeSize = 5)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
print(resultFisher)

count = 24
goTermTab<-GenTable(GOdata, classicFisher = resultFisher, topNodes = count )

write.table(goTermTab, file = outputFileDF, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

