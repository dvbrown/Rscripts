library(affy)
library(limma)
library(survival)
setwd('~/Documents/CREB/gliomaMicroarray/Gravendeel/analysis/')
source('~/Documents/Rscripts/120524-annotate-probes.R')
rawData = ReadAffy()

subtract = function(x,y) {
  z = x - y
  return(z)
}
divide = function(x,y) {
  z = x / y
  return(z)
}
zTransform = function(matrix, rowMean, rowStdDev) { #convert to the unitless z-score based on a normal distribution
    meanT = apply(matrix, 2, FUN=subtract, rowMean)
    sdT = apply(meanT, 2, FUN=divide, rowStdDev)
    return (sdT)
}

extractDataClinical = function(clinicalData) { #A function to get useful clinical data and match rownames to expression
  sur = clinicalData[,,]
  barcode = as.character(clinicalData[,1])
  barcode = gsub('$', '.CEL', barcode) #append .CEL to the first row to match the array filenames
  sur = cbind(barcode, sur)
  row.names(sur) = sur$barcode
  return (sur)
}
betterGeneScore = function(GeneExpressionMatrix, upRegGenes, downRegGenes) {
  #takes as input a z-score transformed gene expression matrix, and 2 dataFrames of upreg and downreg signature genes
  #subset the gene expression matrix with the signature
  upGenes = GeneExpressionMatrix[names(upRegGenes),]
  downGenes = GeneExpressionMatrix[names(downRegGenes),]
  upSig = colMeans(upGenes, na.rm=T)
  downSig = colMeans(downGenes, na.rm=T)
  score = upSig - downSig
  return (score)
}
bindSignatureToSurvival = function(geneScore, clinical) {
  #bind the computed gene signature score to clinical data for analysis
  geneSet = as.data.frame(geneScore)
  set = intersect(row.names(geneSet), row.names(clinical)) #get the overlap
  data = cbind(geneSet[set,], clinical[set,]) #bind the 2 datasets together based on overlap
  colnames(data) = c('sigScore', 'barcode', 'filename','number', 'sex', 'grade','status','survivalMonth', 'survivalYear')
  row.names(data) = data$barcode
  data = data[,c(1,2,4,5,6,7,8,9)]
  data$survivalMonth = as.numeric(data$survivalMonth) 
  data$survivalYear = as.numeric(data$survivalYear)
  return(data)
}
boxData = exprs(rawData[,100:110])
boxplot(boxData, col=rainbow(11), main='Gravendeel unNormalised data')

data = rma(rawData, verbose=T)
boxNorm = exprs(data[,100:110])
boxplot(boxNorm, main='Gravendeel RMAnormalised data', col=rainbow(11))

exprData = as.data.frame(exprs(data))
exprData = merge.data.frame(genelist, exprData, by.x='GeneID', by.y='row.names')
write.table(exprData, './analysis/130121_duplicateGenes.txt', sep='\t', row.names=F)
#have to run the filterProbes.py script to remove duplicate entries, keeping only the probe with highest sum. Remember to put the headers on top! 

#get the gene expression matrix in the right format
gse = read.delim('130122_filterGenes.txt', sep='\t')
row.names(gse) = gse$GeneSymbol
gse = gse[,2:233]
gem = as.matrix(gse)
rowMean = rowMeans(gem)
rowStdDev = apply(gem, 1, sd)
#compute the z-scores for the dataFrame
zScore = zTransform(gem, rowMean, rowStdDev)

write.table(zScore, '130123_zScoreGravendeelTransformeded.txt', sep='\t')

#CREB signature change this!
signature = read.delim('~/Documents/CREB/ChIPseqENCODE/doc/130526_ENCODEsignature.txt')
signature = as.character(signature[,1]) #the score column ranks the number of sites in the given promoter
signatureOut = '~/Documents/CREB/CREBsignatures/GravandeelBoundData/130130_LembergerCREBkoOnlyTop50.txt'

#format the clinical data from Gravendeel
clinical = read.delim('~/Documents/CREB/gliomaMicroarray/Gravendeel/analysis/design.txt')
clinical = extractDataClinical(clinical)

gravendeel = as.data.frame(zScore)
subsetGrav = gravendeel[signature,] #subset the data with the gene list. This is the signature     
upDown = (sort(apply(subsetGrav, 1, median))) #get the median expression score to define up and down genes
#subset the gene score for upregulated and downregulated genes
up = upDown[upDown > 0] #upregulaed genes input for the geneSignature score function
down = upDown[upDown < 0]

#compute the geneScore
geneScore = betterGeneScore(gravendeel, up, down)
#attach clinical data to the signature score
data = bindSignatureToSurvival(geneScore, clinical)

#generate a column listing above the kth percentile
percentile = quantile(data$sigScore, probs=0.50, names=T)
data$percentile = ifelse(data$sigScore >= percentile, 'high', 'low')
write.table(data, signatureOut, sep='\t')

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(data$survivalMonth, event=data$status)
sur.fit =survfit(data.surv~data$percentile)
plot(sur.fit, main='CREB ChIP Gravandeel',ylab='Survival probability',xlab='survival(months)', col=c('red','blue'),xlim=c(0,750))
legend('topright', c('CREB score > 50th, n=102', 'CREB score < 50th, n=94'), col=c('red', 'blue'),lwd=1, cex=0.6)
summary(data.surv)
#test for a difference between curves
test = survdiff(data.surv~data$percentile)
test

text(locator(1),labels='p=0.0136', cex=1) #add the p-value to the graph

#plot the CREB score for each tumour grade. Go into initaialise R to fix cex, las, mar(9,7,5,2)) and mgp
boxplot(data$sigScore~data$grade, col=rainbow(8), main='CREB target genes Gravendeel glioma data set',xlab='Tumour Grade',
        ylab='CREB signature score', par(cex=1.25,las=2))

library(RColorBrewer)
matGrav = as.matrix(subsetGrav)
heatmap(matGrav, margins=c(7,5),cexRow=0.5, Colv=data$grade, labCol=NA, Rowv=NA,xlab='Patients', col=brewer.pal(9,"YlOrRd"),
        ylab='CREB target gene set', main='Gravandeel glioma dataset CREB signature')