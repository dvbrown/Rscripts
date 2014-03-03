library(DriverNet)

setwd('~/Documents/public-datasets/RNA-seq/driverNet/')
# Import data
mut = data(samplePatientMutationMatrix)
head(mut)
gene = data(samplePatientOutlierMatrix)
head(gene)

path = data(sampleInfluenceGraph)
head(path)
# This is a two-dimensional binary matrix. The row and column names are gene names.
# If two genes i and j are in the same pathway, sampleInfluenceGraph[i, j] = 1.

geneNames = data(sampleGeneNames)

# The main function to compute drivers
driversList = computeDrivers(samplePatientMutationMatrix,
    samplePatientOutlierMatrix,sampleInfluenceGraph,
    outputFolder=NULL, printToConsole=T)
drivers(driversList)[1:10]

# random permute the gene labels to compute p-values
randomDriversResult = computeRandomizedResult(
    patMutMatrix=samplePatientMutationMatrix,
    patOutMatrix=samplePatientOutlierMatrix,
    influenceGraph=sampleInfluenceGraph,
    geneNameList= sampleGeneNames, outputFolder=NULL,
    printToConsole=FALSE,numberOfRandomTests=20, weight=FALSE,
    purturbGraph=FALSE, purturbData=TRUE)

# Summarize the results
res = resultSummary(driversList, randomDriversResult,
    samplePatientMutationMatrix,sampleInfluenceGraph,
    outputFolder="", printToConsole=FALSE)
res[1:2,]

# The actual events is poorly documented but it appears to list for each driver gene which patient and gene is 
# Misregulated in a list like format
events = actualEvents(driversList)