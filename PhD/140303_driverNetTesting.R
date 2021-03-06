library(DriverNet)
library(reshape)

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

############################################################## Make the mutational matrix sparse like the TruSeq cancer amplicon ############################################################## 

panel = c('ABL1', 'EGFR', 'GNAS', 'MLH1', 'RET', 'AKT1', 'ERBB2', 'HNF1A', 'MPL', 'SMAD4', 'ALK', 'ERBB4', 'HRAS', 'NOTCH1', 'SMARCB1', 'APC', 'FBXW7', 'IDH1', 'NPM1',
          'SMO', 'ATM', 'FGFR1', 'JAK2', 'NRAS', 'SRC', 'BRAF', 'FGFR2', 'JAK3', 'PDGFRA', 'STK11', 'CDH1', 'FGFR3', 'KDR', 'PIK3CA', 'TP53', 'CDKN2A', 
          'FLT3', 'KIT', 'PTEN', 'VHL', 'CSF1R', 'GNA11', 'KRAS', 'PTPN11', 'CTNNB1', 'GNAQ', 'MET', 'RB1')	

panelMut = samplePatientMutationMatrix[,colnames(samplePatientMutationMatrix) %in% panel]

# The main function to compute drivers
driversList = computeDrivers(panelMut,
                             samplePatientOutlierMatrix,sampleInfluenceGraph,
                             outputFolder=NULL, printToConsole=T)
drivers(driversList)[1:10]

# random permute the gene labels to compute p-values
randomDriversResult = computeRandomizedResult(
    patMutMatrix=panelMut,
    patOutMatrix=samplePatientOutlierMatrix,
    influenceGraph=sampleInfluenceGraph,
    geneNameList= sampleGeneNames, outputFolder=NULL,
    printToConsole=FALSE,numberOfRandomTests=20, weight=FALSE,
    purturbGraph=FALSE, purturbData=TRUE)

# Summarize the results
results = resultSummary(driversList, randomDriversResult,
                    panelMut,sampleInfluenceGraph,
                    outputFolder="", printToConsole=FALSE)
results[1:2,]

# The actual events is poorly documented but it appears to list for each driver gene which patient and gene is 
# Misregulated in a list like format
events = actualEvents(driversList)

res
results
# Compared both the results from exome and the targeted panel and conclude that drverNet works alright.
# Don't capture the full extent of events but the p values and other metrics are the same.