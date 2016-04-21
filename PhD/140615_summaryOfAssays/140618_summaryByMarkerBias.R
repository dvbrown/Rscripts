setwd('~/Documents/Cell_biology/proliferation/Resazurin/140614_summary/analysisByMarkerBias/')
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

calcDMSOcontrol = function(dataFrame) {
    vehicle = dataFrame[dataFrame$treatment %in% 'DMSO',]
    tmz = dataFrame[dataFrame$treatment %in% 'TMZ',]  
    tmz$rep1 = tmz$rep1 / vehicle$rep1
    tmz$rep2 = tmz$rep2 / vehicle$rep2
    tmz$rep3 = tmz$rep3 / vehicle$rep3
    tmz$mean = rowMeans(tmz[,c(4:6)], na.rm=T)
    tmz$sd = apply(tmz[,c(4:6)], 1, sd, na.rm=T)
    return (tmz)
}

calcProlifNormalised = function(dataFrame) {
    negative = dataFrame[dataFrame$cd133 %in% 'CD133_neg',]
    positive = dataFrame[dataFrame$cd133 %in% 'CD133_pos',]
    positive$rep1 = positive$rep1 / negative$rep1
    positive$rep2 = positive$rep2 / negative$rep2
    positive$rep3 = positive$rep3 / negative$rep3
    positive$mean = rowMeans(positive[,c(4:6)], na.rm=T)
    positive$sd = apply(positive[,c(4:6)], 1, sd, na.rm=T)
    return (positive)
}

extractPosNegReplicates = function(dataFrame) {
    neg = dataFrame[dataFrame$cd133 %in% 'neg',c(1,2,3,9)]
    pos = dataFrame[dataFrame$cd133 %in% 'pos',c(1,2,3,9)]
    negMean = mean(neg$dmsoCorrected)
    negSD = sd(neg$dmsoCorrected)
    posMean = mean(pos$dmsoCorrected)
    posSD = sd(pos$dmsoCorrected)
    negSummary = c(negMean, negSD)
    posSummary = c(posMean, posSD)
    result = rbind(negSummary, posSummary)
    result = as.data.frame(result)
    origin = c('negative', 'positive')
    result = cbind(origin, result)
    colnames(result) = c('origin', 'mean', 'sd')
    return (result)
}

normaliseCD133 <- function (dataFrame) {
    cd133Neg = dataFrame[dataFrame$cd133status %in% 'CD133_neg',]
    cd133Pos = dataFrame[dataFrame$cd133status %in% 'CD133_pos',]
    cd133NegAv = mean(cd133Neg$mean)
    cd133NegSd = sd(cd133Neg$mean) / sqrt(length(cd133Neg$mean))
    cd133PosAv = mean(cd133Pos$mean)
    cd133PosSd = sd(cd133Pos$mean) / sqrt(length(cd133Pos$mean))
    cd133 = as.data.frame(rbind(c(cd133NegAv, cd133NegSd), c(cd133PosAv, cd133PosSd)))
    cd133$cd133 = c('negative', 'positive')
    return (cd133)
}

################################ IO and subsetting ############################################
resazurin = read.delim('~/Documents/Cell_biology/proliferation/Resazurin/140614_summary/140614_day7meanSD.txt')
invasion = read.delim('~/Documents/Cell_biology/microscopy/invasion/140615_summary/140617_summary.rep.txt')
clinical = read.delim('~/Documents/Cell_biology/140618_cloneProgressAssays.txt')
invasion$clone = c('004','004','004','004','020','020','020','034','034','034',
                '034','035','035','035','035','039','039','039','039','041','041','041', '039')
invasion$clone = as.factor(invasion$clone)

# I subset this data for which the dominat population was in the FACS data
dominantGrowth = read.delim('~/Documents/Cell_biology/proliferation/Resazurin//140614_summary/analysisByMarkerBias/140618_growthSubsetDominantSubtype.txt', 
                            colClasses=c('factor','factor', 'factor', 'numeric','numeric', 'numeric', 'numeric', 'numeric', 'factor'))
dominantGrowth$sample = paste(dominantGrowth$clone, dominantGrowth$treatment)

# Add the population bias to the data
resazurinC = merge.data.frame(resazurin, clinical, by.x='clone', by.y='clone')
resazurinC$sample = paste(resazurinC$clone, resazurinC$cd133status)
invasionC = merge.data.frame(invasion, clinical, by.x='clone', by.y='clone')
dominantGrowthC = merge.data.frame(dominantGrowth, clinical, by.x='clone', by.y='clone')

# Subset the data for what stain I used
dataCD133 = resazurinC[resazurinC$cd133status %in% c("CD133_neg", "CD133_pos"),]
dataMatched = dataCD133[!dataCD133$clone %in% c("030a", "034a"),]
dataDoubleStain = resazurinC[!resazurinC$cd133status %in% c("CD133_neg", "CD133_pos"),]

################################  Plot the clinical data by receptor bias ################################ 
clinicalP = ggplot(data=clinical[!is.na(clinical$populationBias),], 
                   aes(x=clone, y=Survival_mo, fill=populationBias)) + 
    scale_fill_manual(values=c("yellow", "skyblue3", "firebrick")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    xlab("Patient") + ylab("Survival (months)") +
    ggtitle("Survival by CD44/ CD133 phenotyping") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


clinicalP2 = ggplot(data=clinical[!is.na(clinical$populationBias),], 
                   aes(x=clone, y=Age, fill=populationBias)) + 
    scale_fill_manual(values=c("yellow", "skyblue3", "firebrick")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    xlab("Patient") + ylab("Age (years)") +
    ggtitle("Age of GBM presentation by CD44/ CD133 phenotyping") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(clinicalP, clinicalP2)

###################################### Regression to CD44 and CD133 #################################################
par(mfrow=c(2,1))
plot(clinical$CD44_percent, clinical$Survival_mo, main="CD44 correlation with survival" , xlab="CD44 expression of GPSC from patient", ylab="Survival (months)")
abline(lm(Survival_mo ~ CD44_percent, clinical), col= "red")
plot(clinical$CD133_percent, clinical$Survival_mo, main="CD133 correlation with survival", xlab="CD133 expression of GPSC from patient", ylab="Survival (months)")
abline(lm(Survival_mo ~ CD133_percent, clinical), col= "red")
par(mfrow=c(1,1))


###################################### Groth plot by FACS bias #################################################

growthP = ggplot(data=dataMatched[dataMatched$treatment %in% 'DMSO',], 
                     aes(x=sample, y=mean, fill=populationBias.x)) + 
    scale_fill_manual(values=c("yellow", "skyblue3", "firebrick")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Patient") + ylab("Fluorescent intensity") +
    ggtitle("Comparing Growth at day 7 \nby CD133 status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
growthP

growthDom = ggplot(data=dominantGrowthC[dominantGrowthC$treatment %in% 'DMSO',], 
                 aes(x=clone, y=mean, fill=populationBias)) + 
    scale_fill_manual(values=c("yellow", "skyblue3", "firebrick")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Patient") + ylab("Fluorescent intensity") +
    ggtitle("Comparing Growth at day 7 \nby CD133 status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
growthDom

day7TMZ = calcDMSOcontrol(dataMatched)
tmzPlot7 = ggplot(data=day7TMZ, aes(x=sample, y=mean, fill=populationBias.x)) + 
    scale_fill_manual(values=c("yellow", "skyblue3", "firebrick")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Patient") + ylab("Cell number relative to DMSO control") +
    ggtitle("Comparing temozolomide sensitivty at day 7\nby CD133 status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
tmzPlot7

day7TMZDom = calcDMSOcontrol(dominantGrowthC)
tmzPlot7D = ggplot(data=day7TMZDom, aes(x=clone, y=mean, fill=populationBias)) + 
    scale_fill_manual(values=c("yellow", "skyblue3", "firebrick")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Patient") + ylab("Cell number relative to DMSO control") +
    ggtitle("Comparing temozolomide sensitivty at day 7\nby CD133 status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
tmzPlot7D