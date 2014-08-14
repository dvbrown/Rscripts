# This was downloaded from the cancer browser which in turn derives from the cancer genome atlas
library(gplots)
source("~/Documents/Rscripts/120704-sortDataFrame.R")

colorVerhaakSubtype <- function (dataFrame) {
    # Colour in the molecular subtypes for heatmap plotting
    # Takes as input a data frame contatining the molecular subtype
    dataFrame$colours = "black"
    dataFrame$colours[dataFrame$GeneExp_Subtype == "Proneural"] = "purple"
    dataFrame$colours[dataFrame$GeneExp_Subtype == "Neural"] = "green"
    dataFrame$colours[dataFrame$GeneExp_Subtype == "Classical"] = "blue"
    dataFrame$colours[dataFrame$GeneExp_Subtype == "Mesenchymal"] = "red"
    dataFrame$colours[dataFrame$G_CIMP_STATUS == "G-CIMP"] = "pink"
    return (dataFrame)
}

colorMySubtype <- function (dataFrame) {
    # Colour in the molecular subtypes for heatmap plotting
    # Takes as input a data frame contatining the molecular subtype
    dataFrame$colours = "black"
    dataFrame$colours[dataFrame$subtype == "CD133"] = "red"
    dataFrame$colours[dataFrame$subtype == "CD44"] = "blue"
    # dataFrame$colours[dataFrame$G_CIMP_STATUS == "G-CIMP"] = "pink"
    return (dataFrame)
}

setwd('~/Documents/public-datasets/cancerBrowser/TCGA_GBM_mutation-2014-05-02/')
rawData = read.delim('genomicMatrix')
# Start with the RNAseq patitents but I have to write the Agilent patients to file
patientSubtype = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/survivalTables/140606_RNAseq_SurvivalBoundData.txt", row.names=1)
#rawData[1,1] = 'blank'

# Remove cases for which there has been no gene assigned
genes = rawData[,1]
rawData = rawData[,c(2:292)]
cases = complete.cases(genes)
data = rawData[cases,]
row.names(data) = genes[cases]

# Subset the subtyped Patietns that exist for the mutations
matched = intersect(row.names(patientSubtype), colnames(data))
dataSubtype = patientSubtype[matched,]
dataMutation = data[,matched]
dataM = as.matrix(dataMutation)

################################# Take only the top 1% mutated genes, which is 3 mutation anyway ##############################
totalMuts = rowSums(dataM)
cutOffMuts = quantile.default(totalMuts, probs=0.99)

dataPresent = dataM[totalMuts >= cutOffMuts,]
totalMuts = rowSums(dataPresent)
toBsorted = as.data.frame(cbind(dataPresent, totalMuts))
 
# Sort the dataframe according to the hightest number of mutations first
dataSort = sort.dataframe(toBsorted, 150, highFirst=T)[,c(1:149)]
dataP = as.matrix(dataSort)
head(dataP)

# Make a heatmap where the input is a true false mutation matrix
heatmap.2(dataP, cexRow=0.75, main="HeatMap of mutations", 
          Colv=dataSubtype$colours, keysize=1, trace="none", col=c('white', 'black'), density.info="none", dendrogram="both", 
          ColSideColors=as.character(dataSubtype$colours), labRow=row.names(dataPresent), 
          xlab="Samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5))
