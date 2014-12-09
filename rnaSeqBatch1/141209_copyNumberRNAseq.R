# Take the gene expression of RNA seq batch1 and attempt to get copy number
library(biomaRt)
library(annotate)
library(org.Hs.eg.db)

annotateIds = function(geneList)  { 
    # Given a gene list return the ENSEMBL information and ID mappings
    ensembl_genes = row.names(geneList)
    list=getBM(
        filters="ensembl_gene_id", 
        attributes=c("ensembl_gene_id", "external_gene_name", "description", "chromosome_name", "start_position", "end_position"),
        values= ensembl_genes,
        mart= mart)
    return (list)
}

setwd("~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/copyNumber/")
# cpm = read.delim("../GLMedgeR/131021_normalisedCPM.txt", row.names=1)
# head(cpm)
# 
# # Load the genomic coordinates in a bed like format
# mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# 
# bedLike = annotateIds(cpm)
# head(bedLike)
# bedCpm = merge.data.frame(bedLike, cpm, by.x="ensembl_gene_id", by.y=0)
# write.table(bedCpm, "141209_bedLike.txt", sep='\t')

bedLike = read.delim("141209_bedLike.txt")