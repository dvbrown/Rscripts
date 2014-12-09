# Take the gene expression of RNA seq batch1 and attempt to get copy number
library(biomaRt)
library(annotate)
library(org.Hs.eg.db)
library(ggbio)
library(GenomicRanges)

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
# Subset only vlaid chromosome names
chr = as.factor(c(1:22, "X", "Y"))
bedLike = bedLike[bedLike$chromosome_name %in% chr,]
droplevels(bedLike$chromosome_name)
seqNames = as.character(bedLike$chromosome_name)

# Implement the copy number conversion from Patel et al
# Sort by chromosome and start position
bedLike = bedLike[order(bedLike$chromosome, bedLike$start_position),]

moveAv <- function(x,n=100){
    return (filter(x,rep(1/n,n), sides=2))
}
moveAv(bedLike$GIC_011)

head(bedLike)
gr = GRanges(seqnames = Rle(seqNames),
            ranges = IRanges(start=bedLike$start_position, end=bedLike$end_position),
            gem = bedLike[,c(7:12)],
            name = bedLike$external_gene_name)
gr