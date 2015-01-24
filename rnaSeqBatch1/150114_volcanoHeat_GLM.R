library(edgeR)
library(ggplot2)
library(Heatplus)
library(RColorBrewer)
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/140203_facsBatch/')

#finalResult = read.delim('GLMedgeR/131021_shortVSlong.txt')
finalResult = read.delim('140203_shortVSlong.txt')
#sigGenes=read.delim('GLMedgeR/131021_shortVSlongDEgenes.txt')
sigGenes=read.delim('140203_shortVSlongDEgenes.txt')

require(ggplot2)
##Highlight genes that have an absolute fold change > 1 and a p-value < Bonferroni cut-off
finalResult$threshold = as.factor(ifelse(abs(finalResult$logFC) > 1 & finalResult$FDR < 0.1, "red", "blue"))

g = ggplot(data=finalResult, aes(x=logFC, y=-log10(FDR), 
                         label=external_gene_id, size=0.2)) + 
        geom_point(shape=19, alpha=0.5, size=1.25, color=finalResult$threshold) +
        xlab("log2 fold change") + ylab("-log10 FDR adjusted p-value") +
        ggtitle("Differential gene expression short vs long term survivors\nPDGCs n=6")  + # Set title
        theme_bw(base_size=18)

#subset gene names for only significant genes
dd_text = finalResult[(abs(finalResult$logFC) > 1) & (finalResult$FDR < 0.00001),]
#add text to volcano
g + geom_text(data = dd_text, size=5, aes(x=logFC, y=-log10(FDR),
                                  label=external_gene_id), colour="black") +
  theme(axis.title.x = element_text(size=20)) + 
  theme(axis.title.y = element_text(size=20)) +
  theme(title=element_text(size=20)) 