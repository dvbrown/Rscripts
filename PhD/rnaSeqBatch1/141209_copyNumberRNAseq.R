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

moveAv <- function(x,n=350){
    # X is the time series a numeric vector
    # n is the window size to use. 100 is from Patel et al 2014 single cell GBM
    # Multily by 8 since I HAVE 8 times as many genes
    # The circular argumenat doesn't appear to make a difference
    return (filter(x,rep(1/n,n), sides=2, circular=T))
}

setwd("~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/copyNumber/")

####################### Generate a bed like file ############################
cpm = read.delim("../GLMedgeR/131021_normalisedCPM.txt", row.names=1)
head(cpm)

# Load the genomic coordinates in a bed like format
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

bedLike = annotateIds(cpm)
head(bedLike)
bedCpm = merge.data.frame(bedLike, cpm, by.x="ensembl_gene_id", by.y=0)
write.table(bedCpm, "141212_bedLike.txt", sep='\t')

#################### Center the expression levels and take the moving average ######################
bedLike = read.delim("141212_bedLike.txt")
# Subset only vlaid chromosome names
chr = as.factor(c(1:22, "X", "Y"))
bedLike = bedLike[bedLike$chromosome_name %in% chr,]
droplevels(bedLike$chromosome_name)

# Implement the copy number conversion from Patel et al
# Sort by chromosome and start position
bedLike = bedLike[order(bedLike$chromosome, bedLike$start_position),]
copyNum = bedLike[,]
copyNum$GIC_011 = as.numeric(moveAv(copyNum$GIC_011))
copyNum$GIC_020 =  as.numeric(moveAv(copyNum$GIC_020))
copyNum$GIC_034 =  as.numeric(moveAv(copyNum$GIC_034))
copyNum$GIC_035 =  as.numeric(moveAv(copyNum$GIC_035))
copyNum$GIC_039 =  as.numeric(moveAv(copyNum$GIC_039))
copyNum$GIC_041 =  as.numeric(moveAv(copyNum$GIC_041))
copyNum$begining = paste(copyNum$chromosome_name, copyNum$start_position, sep="_")

# Extract and z center the copy number
copMat = as.matrix(copyNum[,c(7:12)])
copMat = scale(copMat)
# Remove NA columns
copyNum[,c(7:12)] = as.data.frame(copMat)
# copyNum = copyNum[c(400:48112),] This is where circular = FALSE

write.table(copyNum, "150116_bedLike_average800.txt", sep='\t')

# 150118 the useful part of the script ends here.