# This was downloaded from the cancer browser which in turn derives from the cancer genome atlas
library(gtools)

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
patientMutation = patientSubtype[matched,]