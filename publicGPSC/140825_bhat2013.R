library(limma)
library(affy)

setwd('~/Documents/public-datasets/GPSC_subgroups/bhat2013/GSE49009_RAW/')
list.files()

rawData = ReadAffy()
dm = readTargets('../designMatrix.txt')

par(mfrow=(c(2,1)))
# boxData = exprs(rawData)
colnames(boxData) = dm$cellLine
boxplot(rawData, col=rainbow(21), main='Bhat unNormalised data', las=1, mar=c(5,5,5,2))
# data = rma(rawData, verbose=T)
# norm = exprs(data)
colnames(norm) = dm$cellLine
boxplot(norm, col=rainbow(21), main='Bhat normalised data', las=2, mar=c(5,5,5,2))
par(mfrow=(c(1,1)))

pca = princomp(norm)
pcaDf = as.data.frame(cbind(pca$loadings[,1], pca$loadings[,2]))
pcaDf = cbind(pcaDf, dm)[c(1:17),]
dd_text = dm$cellLine

g = ggplot(data=pcaDf, aes(x=V1, y=V2, color=Subtype)) + 
    geom_point(shape=19) + #geom_smooth(method=lm, colour='red') +
    xlab("PC1") + ylab("PC2") + # Set axis labels
    ggtitle("Principal component analysis Bhat 2013") +  # Set title
    theme_bw(base_size=18)

g + geom_text(data = pcaDf, aes(x=V1, y=V2,
        label=dd_text[1:17], size=0.2), colour="black")

g2 = ggplot(data=pcaDf, aes(x=V1, y=V2, color=Source)) + 
    geom_point(shape=19) + #geom_smooth(method=lm, colour='red') +
    xlab("PC1") + ylab("PC2") + # Set axis labels
    ggtitle("Principal component analysis Bhat 2013") +  # Set title
    theme_bw(base_size=18)

g2 + geom_text(data = pcaDf, aes(x=V1, y=V2,
                                label=dd_text[1:17], size=0.2), colour="black")