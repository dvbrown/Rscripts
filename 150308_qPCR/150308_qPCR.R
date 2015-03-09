library(reshape)
library(plyr)
library(ggplot2)
source('./multiplot.R')
source('../qPCRFunctions.R')

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

function(geneOfInterest, sampleOfInterest='020_N', houseKeepingGene='GAPDH', referenceSample='020_N', data=rawData) {
    sampleHouse = paste(sampleOfInterest, houseKeepingGene)
    SampleGene = paste(sampleOfInterest, geneOfInterest)
    # Extract the Cp of the hosue keeping gene
    houseCp = data[sampleHouse, 'meanCP']
    geneCp = data[SampleGene, 'meanCP']
    # dCt calculation for the sample of interest
    dCt = houseCp - geneCp
    # Extract the meanCP for the reference sample. First get the index of the housekeeping gene, then the gene of interest
    refDctRowHouse = paste(referenceSample, houseKeepingGene)
    refDctRowGene = paste(referenceSample, geneOfInterest)
    # Calculate dCt for the reference sample
    referenceSample_dCt = data[refDctRowHouse, 'meanCP'] - data[refDctRowGene, 'meanCP']
    # Calculate ddCt
    ddCtNotSquared = dCt - referenceSample_dCt
    ddCt = 2^ddCtNotSquared
    return (ddCt)
}

c035 = datWide[datWide$PDGC %in% "MU035",]
c035$ddCt = ddCTcalculate(geneOfInterest=c035$Gene, sampleOfInterest=c035$cDNA,
                          houseKeepingGene='GAPDH', referenceSample='MU035_DN', data=c035) 