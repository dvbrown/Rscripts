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
dat$Gene = gsub("POU5F1", "OCT4", dat$Gene)
dat$Gene = gsub("PROM1", "CD133", dat$Gene)

# Take the mean of duplicated measurements
datWide = ddply(dat, .(Sample, Gene, PDGC, Subpopulation), summarise, meanCp = mean(MeanCP, na.rm=T))

# Calculate ddCT

getddCt = function(dataFrame, sampleInt="MU035_CD133_neg") 
{
    # Extract the genes that were measured by both samples.
    samp = dataFrame[dataFrame$Sample %in% sampleInt,"Gene"]
    cd133neg = gsub("_CD133_pos", "_CD133_neg", sampleInt)
    ref = dataFrame[dataFrame$Sample %in% cd133neg, "Gene"]
    genes = intersect(samp, ref)
    
    # Subset the shared genes and sort
    subDat = dataFrame[dataFrame$Sample %in% c(sampleInt, cd133neg),]
    subDat = subDat[subDat$Gene %in% genes,]
    subDat = subDat[order(subDat$Sample, subDat$Gene),]
    row.names(subDat) = paste(subDat$Sample, subDat$Gene)
    
    # Call the ddCt functions
    subDat$ddCT = ddCTcalculate(geneOfInterest=subDat$Gene, sampleOfInterest=subDat$Sample,
                                houseKeepingGene='GAPDH', referenceSample=cd133neg, data=subDat)
    subDat$foldChange = foldChangecalculate(geneOfInterest=subDat$Gene, sampleOfInterest=subDat$Sample,
                                            houseKeepingGene='GAPDH', referenceSample=cd133neg, data=subDat)
    
    # Remove the reference samples which will be 0 anyway.
    result = subDat[subDat$Sample %in% sampleInt,]
    return (subDat)
    #return (result)
}

############   Calculate ddCT ############  
dataFrame = datWide
sampleInt = "MU035_CD133_pos"

mu035_P = getddCt(datWide, "MU035_CD133_pos")
mu011_P = getddCt(datWide, "MU011_CD133_pos")
mu020_P = getddCt(datWide, "MU020_CD133_pos")
mu030_P = getddCt(datWide, "MU030_CD133_pos")
mu030a_P = getddCt(datWide, "MU030a_CD133_pos")
mu041_P = getddCt(datWide, "MU041_CD133_pos")

ddCt = rbind(mu011_P, mu020_P, mu030_P, mu030a_P, mu035_P, mu041_P)
# write.csv(ddCt,"150317_ddCtResults.csv", row.names=F)

rm(mu011_P, mu020_P, mu030_P, mu030a_P, mu035_P, mu041_P)
# Get rid of the genes only measured on ES cells
ddCt = ddCt[!ddCt$Gene %in% c("ATP5G3", "B2M", "CREB1",  "GRIA2", #"GAPDH",
                              "HAPLN1", "MGMT", "REST"),]

orderedDat = ddCt[order(ddCt$Gene, ddCt$PDGC),]

# Make a boxplot of ddCT valuse
ggplot(data=ddCt, aes(x=Gene, y=ddCT)) +
    geom_boxplot() + geom_point(aes(colour=Sample), size =3, alpha=0.7,  position = position_jitter(w = 0.2)) +
    ggtitle("qPCR Summary") + geom_hline(yintercept=0, colour="red") +
    xlab("Gene") + ylab("ddCt relative to CD133 negative") +
    scale_y_continuous(breaks = round(seq(-6, 6, by = 1),1)) +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

######## Summarise by biological replicates ######## 

bioRep = ddply(ddCt, .(Subpopulation, Gene), summarise, meanddCt = mean(ddCT, na.rm=T), 
               sdddCt = sd(ddCT, na.rm=T), reps=length(ddCT))
bioRep$seddCt = bioRep$sdddCt / (bioRep$reps)
write.csv(bioRep, "150321_biologicalRepCD133.csv")

# Plot the biological replicates
ggplot(data=bioRep, aes(x=Gene, y=meanddCt, fill=Subpopulation)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    ggtitle("qPCR biological Replicates (n = 5 - 6)") +  scale_fill_manual(values=c("blue", "gold1")) + 
    geom_errorbar(aes(ymin=meanddCt-seddCt, ymax=meanddCt+seddCt), width=.2, position=position_dodge(0.9)) +
    xlab("Gene") + ylab("ddCT relative to H9 NSC") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))