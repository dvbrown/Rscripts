library(reshape)
library(plyr)
library(ggplot2)
#source('~/Documents/Rscripts/multiplot.R')

ddCTcalculate = function(geneOfInterest, sampleOfInterest='020_N', houseKeepingGene='GAPDH', referenceSample='020_N', data=rawData) 
{
    sampleHouse = paste(sampleOfInterest, houseKeepingGene)
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
    sampleHouse = paste(sampleOfInterest, houseKeepingGene)
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

plotPDGC = function(dataFrame, pdgc="MU035") {
    # This function takes a dataFrane of some data and returns a ggplot object specific to the cell line supplied as a string
    extract = dataFrame[,"PDGC"] %in% pdgc
    datPlot = dataFrame[extract,]
    # The colour of the plots
    cols = c("orange", "forestgreen", "darkblue", "darkred")
    p = ggplot(data=datPlot, aes(x=Gene, y=ddCT, fill=Subpopulation)) +
        geom_bar(stat="identity", position=position_dodge(), colour="black") + 
        ggtitle(paste("qPCR", pdgc, sep=" ")) +  scale_fill_manual(values=cols) + 
        xlab("Gene") + ylab("ddCt") +
        theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    return (p) }

setwd("150316_qPCRcd133Summary/")

dat = read.csv("dat/150316_summaryEdit.csv", row.names=1)
dat$Sample = paste(dat$PDGC, dat$Subpopulation, sep="_")

# Take the mean of duplicated measurements
datWide = ddply(dat, .(Sample, Gene, PDGC, Subpopulation), summarise, meanCp = mean(MeanCP, na.rm=T))

# Calculate ddCT

getddCt = function(dataFrame, sampleInt="MU035_CD133_neg") 
    {
    # Extract the genes that were measured by both samples.
    samp = droplevels(dataFrame[dataFrame$Sample %in% sampleInt,"Gene"])
    ref = droplevels(dataFrame[dataFrame$Sample %in% 'H9_NSC_NA', "Gene"])
    genes = intersect(samp, ref)
    
    # Subset the shared genes and sort
    subDat = dataFrame[dataFrame$Sample %in% c(sampleInt, "H9_NSC_NA"),]
    subDat = subDat[subDat$Gene %in% genes,]
    subDat = subDat[order(subDat$Sample, subDat$Gene),]
    row.names(subDat) = paste(subDat$Sample, subDat$Gene)
    
    # Call the ddCt functions
    subDat$ddCT = ddCTcalculate(geneOfInterest=subDat$Gene, sampleOfInterest=subDat$Sample,
                                houseKeepingGene='GAPDH', referenceSample='H9_NSC_NA', data=subDat)
    subDat$foldChange = foldChangecalculate(geneOfInterest=subDat$Gene, sampleOfInterest=subDat$Sample,
                                houseKeepingGene='GAPDH', referenceSample='H9_NSC_NA', data=subDat)
    
    # Remove the reference samples which will be 0 anyway.
    result = subDat[subDat$Sample %in% sampleInt,]
    #return (subDat)
    return (result)
}

############   Calculate ddCT ############  
mu035_N = getddCt(datWide, "MU035_CD133_neg")
mu035_P = getddCt(datWide, "MU035_CD133_pos")
mu011_N = getddCt(datWide, "MU011_CD133_neg")
mu011_P = getddCt(datWide, "MU011_CD133_pos")
mu020_N = getddCt(datWide, "MU020_CD133_neg")
mu020_P = getddCt(datWide, "MU020_CD133_pos")
mu030_N = getddCt(datWide, "MU030_CD133_neg")
mu030_P = getddCt(datWide, "MU030_CD133_pos")
mu030a_N = getddCt(datWide, "MU030a_CD133_neg")
mu030a_P = getddCt(datWide, "MU030a_CD133_pos")
mu041_N = getddCt(datWide, "MU041_CD133_neg")
mu041_P = getddCt(datWide, "MU041_CD133_pos")
h9ES = getddCt(datWide, "H9_ES_NA")

ddCt = rbind(mu011_N, mu011_P, mu020_N, mu020_P,
             mu030_N, mu030_P, mu030a_N, mu030a_P,
             mu035_N, mu035_P, mu041_N, mu041_P, h9ES)
write.csv(ddCt,"150317_ddCtResults.csv", row.names=F)

rm(mu011_N, mu011_P, mu020_N, mu020_P, mu030_N, mu030_P, h9ES,
   mu030a_N, mu030a_P, mu035_N, mu035_P, mu041_N, mu041_P)

ggplot(data=ddCt, aes(x=Gene, y=ddCT, fill=Sample)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    ggtitle("qPCR individual samples") +  #scale_fill_manual(values=cols) + 
    xlab("Gene") + ylab("ddCt relative to H9 NSC") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

######## Summarise by bioloical replicates ######## 

