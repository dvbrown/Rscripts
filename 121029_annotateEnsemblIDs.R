library(biomaRt)

mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

annotateIds = function(geneList)  { #Given a gene list return the ENSEMBL information and ID mappings
  ensembl_genes = geneList$id
  list=getBM(
           filters= "ensembl_gene_id", 
           attributes= c("ensembl_gene_id", "external_gene_id", "entrezgene", "description"),
           values= ensembl_genes,
           mart= mart)
  return (list)
}

#retrieve gene names map to DE test
annResult = annotateIds(result)
finalResult = merge.data.frame(result, annResult, by.x='id', by.y='ensembl_gene_id')

annResult.unsort.negative = annotateIds(result.unsort.negative)
final.result.unsort.negative =  merge.data.frame(result.unsort.negative, annResult.unsort.negative,by.x='id', by.y='ensembl_gene_id')
  
annResult.unsort.positive = annotateIds(result.unsort.positive)
final.result.unsort.positive = merge.data.frame(result.unsort.positive, annResult.unsort.positive, by.x='id', by.y='ensembl_gene_id')

write.table(finalResult, './121030_correctSamples/121030_CD133lowVSCD133hi.txt', sep='\t')
write.table(final.result.unsort.negative, './121030_correctSamples/121030_unsortVSCD133low.txt', sep='\t')
write.table(final.result.unsort.positive, './121030_correctSamples/121030_unsortVSCD133high.txt', sep='\t')