library(ggplot2)
library(reshape)
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/140203_facsBatch/')

dat = read.delim("140213_normalisedLog_CPM_facs.txt")
genes = read.delim('140203_shortVSlong.txt')[,c(1,2)]

# Join on gene names
datJoin = merge.data.frame(genes, dat, by.x="ensembl_gene_id", by.y=0)

# List the genes which I want to plot expression scores for
gen = c("TUBB3", "PROM1", "GFAP", "NANOG", "NES", "OLIG2", "SOX2", "POU5F1", "LAMB1")

# Subset
datSub = datJoin[datJoin$external_gene_id %in% gen,]
datSub = datSub[c(1:5,7,8),]
row.names(datSub) = datSub$external_gene_id
datSub = datSub[,c(2:8)]

# Melt into for for ggplot
datMelt = melt(datSub, id.vars="external_gene_id")
colnames(datMelt) = c("Gene", "PDGC", "Log_CPM")

# Colour by survival status
datMelt$Status = ""
datMelt[datMelt$PDGC %in% c("GIC_011", "GIC_020", "GIC_034"),]$Status = "Short-term"
datMelt[!datMelt$PDGC %in% c("GIC_011", "GIC_020", "GIC_034"),]$Status = "Long-term"

# Make a boxplot of ddCT valuse
ggplot(data=datMelt, aes(x=Gene, y=Log_CPM)) +
    geom_boxplot() + geom_point(aes(colour=Status), size =3, alpha=0.7,  position = position_jitter(w = 0.2)) +
    scale_colour_manual(values=c("blue", "forestgreen")) + 
    xlab("Gene") + ylab("Log CPM gene expression") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    theme(text = element_text(size=20))