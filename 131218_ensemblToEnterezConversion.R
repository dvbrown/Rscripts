# A helper function to generate a dataframe mapping enterez to ensembl IDs.
# Includes deduplication function

library(annotate)
library(org.Hs.eg.db)

ensembl2enterezConvert = function(genesTestedforDE) {
  ensembl2enterez = select(org.Hs.eg.db, keys=row.names(genesTestedforDE), cols='ENTREZID', keytype='ENSEMBL')
  ensembl2enterez = ensembl2enterez[!is.na(ensembl2enterez),]
  ensembl2enterez = ensembl2enterez[!is.na(ensembl2enterez$ENTREZID),]
  ensembl2enterez = ensembl2enterez[!duplicated(ensembl2enterez$ENTREZID),]
  return (ensembl2enterez)
}