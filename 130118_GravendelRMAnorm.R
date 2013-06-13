library(affy)
library(limma)
library(survival)
setwd('~/Documents/CREB/gliomaMicroarray/Gravendeel/analysis/')
source('~/Documents/Rscripts/120524-annotate-probes.R')
source('~/Documents/Rscripts/121112_subSetTCGAfunctions.R')
rawData = ReadAffy()

zTransform = function(matrixElement, rowMean, rowSD ) { #convert to the unitless z-score based on a normal distribution
  z = (matrixElement - rowMean)/rowSD
  return (z)
}

extractDataClinical = function(clinicalData) { #A function to get useful clinical data and match rownames to expression
  sur = clinicalData[,,]
  barcode = as.character(clinicalData[,1])
  barcode = gsub('$', '.CEL', barcode) #append .CEL to the first row to match the array filenames
  sur = cbind(barcode, sur)
  row.names(sur) = sur$barcode
  return (sur)
}
bindSignatureToSurvival = function(geneScore) {
  #bind the computed gene signature score to clinical data for analysis
  geneSet = as.data.frame(geneScore)
  set = intersect(row.names(geneSet), row.names(clinical)) #get the overlap
  data = cbind(geneSet[set,], clinical[set,]) #bind the 2 datasets together based on overlap
  colnames(data) = c('sigScore', 'barcode', 'filename','number', 'sex', 'grade','status','survivalMonth', 'survivalYear')
  row.names(data) = data$barcode
  data = data[,c(1,2,4,5,6,7,8,9)]
  data$survivalMonth = as.double(data$survivalMonth)
  data$survivalYear = as.double(data$survivalYear)
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
gse$dup = duplicated(gse$GeneSymbol)
gem = subset.data.frame(gse, gse$dup == FALSE)
row.names(gem) = gem$GeneSymbol
gem = gem[,c(2:233)]
gem = as.matrix(gem)

#compute the z-scores for the dataFrame
rowMean = rowMeans(gem)
rowStdDev = apply(gem, 1, sd)
zScore = apply(gem, 2, zTransform, rowMean, rowStdDev) 
write.table(zScore, '130122_zScoreGravendeelTransformed.txt', sep='\t')

#Retrieve the Gravandeel z-score data
gravendeel = read.delim('~/Documents/CREB/gliomaMicroarray/Gravendeel/analysis/130122_zScoreGravendeelTransformed.txt', row.names=1)
gravendeel = as.data.frame(gravendeel)
#CREB signature change this!
signature = read.delim('~/Documents/CREB/PaulList/130121_PI3KMaPKdependantCREBtargetgenes.txt')
signature = as.character(signature[,1]) #the score column ranks the number of sites in the given promoter

#format the clinical data from Gravendeel
clinical = read.delim('~/Documents/CREB/gliomaMicroarray/Gravendeel/analysis/design.txt')
clinical = extractDataClinical(clinical)

subsetGrav = gravendeel[signature,] #subset the data with the gene list. This is the signature     
upDown = (sort(apply(subsetGrav, 1, median))) #get the median expression score to define up and down genes
#subset the gene score for upregulated and downregulated genes
up = upDown[upDown > 0] #upregulaed genes input for the geneSignature score function
down = upDown[upDown < 0]

#compute the geneScore
geneScore = betterGeneScore(gravendeel, up, down)
#attach clinical data to the signature score
data = bindSignatureToSurvival(geneScore)

#generate a column listing above the kth percentile
percentile = quantile(data$sigScore, probs=0.50, names=T)
data$percentile = ifelse(data$sigScore >= percentile, 'high', 'low')

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(data$survivalMonth, event=data$status)
sur.fit =survfit(data.surv~data$percentile)
plot(sur.fit, main='Gravendeel glioma data set',ylab='Survival probability',xlab='survival(months)', col=c('red','blue'),xlim=c(0,250))
legend('topright', c('CREB score > 50th, n=103', 'Stem cell score < 50th, n=93'), col=c('red', 'blue'),lwd=1, cex=0.6)
summary(data.surv)
#test for a difference between curves
test = survdiff(data.surv~data$percentile)
test

text(locator(1),labels='p=0.23', cex=1) #add the p-value to the graph

#plot the CREB score for each tumour grade. Go into initaialise R to fix cex, las, mar(9,7,5,2)) and mgp
boxplot(data$sigScore~data$grade, col=rainbow(8), main='CREB target genes Gravendeel glioma data set',xlab='Tumour Grade',
        ylab='CREB signature score', par(cex=1.25))