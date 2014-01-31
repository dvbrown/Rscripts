setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/limmaResults/')
list.files()
# Read in the files into a list
affyFiles = list.files(pattern='*.affymetrix*')
agilentFiles = list.files(pattern='*.agilent*')
affy = lapply(affyFiles, read.delim, header=T)
names(affy) = affyFiles

agilFiles = lapply(agilentFiles, read.delim, header=T)
names(agilFiles) = agilentFiles

# Build a dataframe of the adjusted p-values from different methods
affyPval = cbind(affy[[1]]$adj.P.Val, affy[[2]]$adj.P.Val, affy[[3]]$adj.P.Val, affy[[4]]$adj.P.Val, affy[[5]]$adj.P.Val)
row.names(affyPval) = affy[[1]]$ID
colnames(affyPval) = affyFiles
affySig = -log10(affyPval)

# Build a dataframe of the adjusted p-values from different methods
agilentPval = cbind(agilFiles[[1]]$adj.P.Val, agilFiles[[2]]$adj.P.Val, agilFiles[[3]]$adj.P.Val)
row.names(agilentPval) = agilFiles[[1]]$ID
colnames(agilentPval) = agilentFiles
agilentSig = -log10(agilentPval)

head(affySig)
head(agilentSig)

# Read in the result of my RNA-seq batch 1
stemCell =  read.delim('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/131021_shortVSlong.txt', row.names=1, stringsAsFactors=F)
genes = !duplicated(stemCell$external_gene_id)
stemCell1 = stemCell[genes,]
row.names(stemCell1) = stemCell1$external_gene_id
stemCellDE = stemCell1[,c(4,5,6,7,8)]

Ag = as.data.frame(agilFiles[[3]])
row.names(Ag) = Ag$ID
Af = as.data.frame(affy[[5]])
row.names(Af) = Af$ID

sigGenes = list(stemCellDE[abs(stemCellDE$logFC) >= 1 & stemCellDE$FDR < 0.1,],
                Ag[abs(Ag$logFC >= 1) & Ag$adj.P.Val < 0.1,],
                Af[abs(Af$logFC >= 1) & Af$adj.P.Val < 0.1,])
names(sigGenes) = c('My RNA-seq', 'Agilent', 'Affymetrix')

# Make some overlap plots
overlap = intersect(Af$ID, stemCell1$external_gene_id)
stemCellcommon = stemCell1[overlap,]
AfCommon = Af[overlap,]

par(mfrow=c(2,1))
plot(AfCommon$logFC, stemCellcommon$logFC, main='logFC stem cell \nvs TCGA bulk GBM',
     ylab='My RNAseq batch 1', xlab='Affymetrix TCGA GBM')
abline(lsfit(AfCommon$logFC, stemCellcommon$logFC), col='red')

plot(-log10(AfCommon$adj.P.Val), -log10(stemCellcommon$FDR), main='Significance stem cell \nvs TCGA bulk GBM',
     ylab='My RNAseq batch 1', xlab='Affymetrix TCGA GBM')
abline(lsfit(-log10(AfCommon$adj.P.Val), -log10(stemCellcommon$FDR)), col='red')


# Make a volcano plot of Affymetrix
Af$threshold = as.factor(abs(Af$logFC > 2 & Af$adj.P.Val < 0.05))
require(ggplot2)

g = ggplot(data=Af, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.80, size=2) +
  theme(legend.position = "none") + ggtitle('Differently expressed genes between short-term survivors\nand > 3 year survivors') +
  #Add some custome limits to the axes here
  #xlim(c(-5, 5)) + ylim(c(0, 50)) +
  xlab(as.character('logFC')) + ylab(as.character('adj.P.Val'))
g

# Add annotation
dd_text = Af[(abs(Af$logFC) > 1) & (Af$adj.P.Val < 0.05),]

# Draw text
g + geom_text(data = dd_text, aes(x=logFC, y=-log10(adj.P.Val),
                                             label=ID, size=0.2), colour="black")