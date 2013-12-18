# A helper function to generate a dataframe mapping enterez to ensembl IDs.
# Includes deduplication function

library(annotate)
library(org.Hs.eg.db)

ensembl2enterezConvert = function(deGenes) {
  ensembl2enterez = select(org.Hs.eg.db, keys=deGenes$ensembl_gene_id, cols='ENTEREZID', keytype='ENSEMBL')
  ensembl2enterez = ensembl2enterez[!is.na(ensembl2enterez),]
  ensembl2enterez = ensembl2enterez[!is.na(ensembl2enterez$ENTEREZID),]
  ensembl2enterez = ensembl2enterez[!duplicated(ensembl2enterez$ENTEREZID),]
  return (ensembl2enterez)
}