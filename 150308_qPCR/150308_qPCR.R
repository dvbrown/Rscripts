library(reshape)
library(plyr)
library(ggplot2)
source('./multiplot.R')
source('../qPCRFunctions.R')

ddCTcalculate = function(geneOfInterest, sampleOfInterest='020_N', houseKeepingGene='GAPDH', referenceSample='020_N', data=rawData) 
{
    sampleHouse = unique(paste(sampleOfInterest, houseKeepingGene))
    sampleGene = paste(sampleOfInterest, geneOfInterest)
    # Extract the Cp of the hosue keeping gene
    houseCp = data[sampleHouse, 'meanCp']
    geneCp = data[sampleGene, 'meanCp']
    # dCt calculation for the sample of interest
    dCt = houseCp - geneCp
    # Extract the meanCP for the reference sample. First get the index of the housekeeping gene, then the gene of interest
    refDctRowHouse = paste(referenceSample, houseKeepingGene)
    refDctRowGene = paste(referenceSample, geneOfInterest)
    # Calculate dCt for the reference sample
    referenceSample_dCt = data[refDctRowHouse, 'meanCp'] - data[refDctRowGene, 'meanCp']
    # Calculate ddCt
    ddCt = dCt - referenceSample_dCt
    return (ddCt)
}

foldChangecalculate = function(geneOfInterest, sampleOfInterest='020_N', houseKeepingGene='GAPDH', referenceSample='020_N', data=rawData) 
{
    sampleHouse = unique(paste(sampleOfInterest, houseKeepingGene))
    sampleGene = paste(sampleOfInterest, geneOfInterest)
    # Extract the Cp of the hosue keeping gene
    houseCp = data[sampleHouse, 'meanCp']
    geneCp = data[sampleGene, 'meanCp']
    # dCt calculation for the sample of interest
    dCt = houseCp - geneCp
    # Extract the meanCP for the reference sample. First get the index of the housekeeping gene, then the gene of interest
    refDctRowHouse = paste(referenceSample, houseKeepingGene)
    refDctRowGene = paste(referenceSample, geneOfInterest)
    # Calculate dCt for the reference sample
    referenceSample_dCt = data[refDctRowHouse, 'meanCp'] - data[refDctRowGene, 'meanCp']
    # Calculate ddCt
    ddCt = dCt - referenceSample_dCt
    foldChange = 2^ddCt
    return (foldChange)
}

setwd("150308_qPCR/")

dat = read.delim("dat/150308_ctValuesMapped.txt", row.names=1)

# Summarise the data by replicate
datWide = ddply(dat, .(cDNA, Gene), summarise, rep1=Cp[1], rep2=Cp[2], meanCp = mean(Cp, na.rm=T), 
                sdCp = sd(Cp), reps=length(Cp))
# Remove the non interesting targets
datWide = datWide[c(1:112),]
# Write to file
write.csv(datWide, "./dat/150308_datWide.csv", row.names=F)

# Plot replicates
correlation = ggplot(data=datWide, aes(x=rep1, y=rep2, color=Gene)) + 
    geom_point(shape=19) + geom_smooth(method=lm, colour='red') +
    xlab("Replicate 1") + ylab("Replicate 2") + # Set axis labels
    ggtitle("Correlation of technical replicates") +  # Set title
    theme_bw(base_size=18)
correlation

#### calculate ddCT ####
# Subset by sample
datWide = cbind(colsplit.factor(datWide$cDNA, split = "_", names=c("PDGC", "Subpopulation")), datWide)
row.names(datWide) = paste(datWide$cDNA, datWide$Gene)

c035 = datWide[datWide$PDGC %in% "MU035",]
c035$ddCT = ddCTcalculate(geneOfInterest=c035$Gene, sampleOfInterest=c035$cDNA,
                          houseKeepingGene='GAPDH', referenceSample='MU035_DN', data=c035)
c035$foldChange = foldChangecalculate(geneOfInterest=c035$Gene, sampleOfInterest=c035$cDNA,
                          houseKeepingGene='GAPDH', referenceSample='MU035_DN', data=c035) 

c041 = datWide[datWide$PDGC %in% "MU041",]
c041$ddCT = ddCTcalculate(geneOfInterest=c041$Gene, sampleOfInterest=c041$cDNA,
                                houseKeepingGene='GAPDH', referenceSample='MU041_DN', data=c041)
c041$foldChange = foldChangecalculate(geneOfInterest=c041$Gene, sampleOfInterest=c041$cDNA,
                                houseKeepingGene='GAPDH', referenceSample='MU041_DN', data=c041) 
# Remove the fcs cells from this analysis
c041 = c041[!c041$Subpopulation %in% "fcs",]

c020 = datWide[datWide$PDGC %in% "MU020",]
c020$ddCT = ddCTcalculate(geneOfInterest=c020$Gene, sampleOfInterest=c020$cDNA,
                          houseKeepingGene='GAPDH', referenceSample='MU020_DN', data=c020)
c020$foldChange = foldChangecalculate(geneOfInterest=c020$Gene, sampleOfInterest=c020$cDNA,
                                      houseKeepingGene='GAPDH', referenceSample='MU020_DN', data=c020) 

markers = datWide[datWide$PDGC %in% c("MU039", "H9", 'MU041'),]
markers$ddCT = ddCTcalculate(geneOfInterest=markers$Gene, sampleOfInterest=markers$cDNA,
                          houseKeepingGene='GAPDH', referenceSample='MU039_neg', data=markers)
markers$foldChange = foldChangecalculate(geneOfInterest=markers$Gene, sampleOfInterest=markers$cDNA,
                                      houseKeepingGene='GAPDH', referenceSample='MU039_neg', data=markers)
# Remove the double sorted samples of MU041 from markers
markers = markers[!is.na(markers$ddCT),]

# Combine the ddCt data
analysed = rbind(c020, c035, c041, markers)

#### Some plots of data ####
# Check the normality of the distribution
hist(analysed$ddCT, breaks="FD")
hist(analysed$foldChange, breaks="FD")