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

moveAv <- function(x,n=100){
    # X is the time series a numeric vector
    # n is the window size to use. 100 is from Patel et al 2014 single cell GBM
    return (filter(x,rep(1/n,n), sides=2, circular=F))
}

setwd("~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/copyNumber/")

####################### Generate a bed like file ############################
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

#################### Center the expression levels and take the moving average ######################
bedLike = read.delim("141209_bedLike.txt")
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
copyNum[,c(7:12)] = as.data.frame(copMat)

write.table(copyNum, "141209_bedLike.txt", sep='\t')

#################### Draw the copy number plot #######################
copyNum = read.delim("141209_bedLike.txt")
copyNum$chromosome_name = paste("chr", copyNum$chromosome_name, sep="")

myLengths =  c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 
                102531392,  90354753,  81195210,  78077248,  59128983,  63025520,  48129895, 
                51304566, 155270560,  59373566)   

seqNames = as.character(copyNum$chromosome_name)
gr = GRanges(seqnames = Rle(seqNames), 
            ranges = IRanges(start=copyNum$start_position, end=copyNum$end_position),
            gem = copyNum[,c(7:12)],
            name = copyNum$external_gene_name,
            startGene=    paste(copyNum$chromosome_name, copyNum$start_position, sep="_"))
gr
seqlengths(gr) = myLengths

smallGR = gr[1:10]
plot<-ggplot(smallGR) + #add basic layer
    geom_point(aes(x=startGene,y=gem.GIC_011)) + #add copy number points to plot
    theme_bw(base_size=24)

autoplot(smallGR)

# and save plot to file
pdf("copy_numbers.png",14,21)
print(plot)
dev.off()