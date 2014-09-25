library(limma)
library(affy)
library(gplots)
library(ggplot2)
source('~/Documents/Rscripts/120704-sortDataFrame.R')
source('~/Documents/Rscripts/multiplot.R')

##################### IO ###########################
setwd('~/Documents/public-datasets/GPSC_subgroups/rath2012_coCultureAstro/GSE37120_RAW/')
list.files()

rawData = ReadAffy()
dm = readTargets('../design.txt')
dm$cellLine = paste(dm[,2], dm$suorce, dm$subpopulation)
myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

par(mfrow=(c(2,1)))
boxData = exprs(rawData)
colnames(boxData) = dm$cellLine
boxplot(rawData, col=rainbow(21), main='Rath unNormalised data', las=1, mar=c(5,5,5,2))
data = rma(rawData, verbose=T)
norm = exprs(data)
colnames(norm) = dm$cellLine
boxplot(norm, col=rainbow(21), main='Rath normalised data', las=2, mar=c(5,5,5,2))
par(mfrow=(c(1,1)))

##################### Remove non-informatives entries ###########################
dm$sample = paste(dm$Source, dm$Subpopulation)
dmCut = dm[c(5,6,11,12,17,18,23,24),]
norm = norm[,c(5,6,11,12,17,18,23,24)]

##################### Render the prinical components ###########################
pca = princomp(norm)
pcaDf = as.data.frame(cbind(pca$loadings[,1], pca$loadings[,2]))
pcaDf = cbind(pcaDf, dmCut)
dd_text = dm$sample

g = ggplot(data=pcaDf, aes(x=V1, y=V2, color=Source)) + 
    geom_point(shape=19) + #geom_smooth(method=lm, colour='red') +
    xlab("PC1") + ylab("PC2") + # Set axis labels
    ggtitle("Principal component analysis Ruth 2012") +  # Set title
    theme_bw(base_size=18)
g
g2 = ggplot(data=pcaDf, aes(x=V1, y=V2, color=Subpopulation)) + 
    geom_point(shape=19) + #geom_smooth(method=lm, colour='red') +
    xlab("PC1") + ylab("PC2") + # Set axis labels
    ggtitle("Principal component analysis Ruth 2012") +  # Set title
    theme_bw(base_size=18)
g2
g3 = ggplot(data=pcaDf, aes(x=V1, y=V2, color=Replicate)) + 
    geom_point(shape=19) + #geom_smooth(method=lm, colour='red') +
    xlab("PC1") + ylab("PC2") + # Set axis labels
    ggtitle("Principal component analysis Ruth 2012") +  # Set title
    theme_bw(base_size=18)
g3
multiplot(g, g2,g3, cols=2)