# Analyse the the clones I was sequencing in batch 1 and compare primary and recurrent
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')
source('~/Documents/Rscripts/qPCRFunctions.R')
################################## IO ################################## 
setwd('~/Documents/RNAdata/qPCRexpt/140206_cd133posNeg/')

# Convert the 384 well plate into the sample label
# Read in qPCR data
plate = transposeLinear('140206_384wellMap.txt')
# Split the transposed file into source and gene
plateMap = splitSampleName(plate)

cp = read.delim('140206_Cps.txt', skip=1)
tm = read.delim('140206_melt.txt', skip=1)

######################################### Filter and prepare data #####################################

data = buildDataFrameForddCT(plateMap, cp)

# Return a datframe with the replicate data alongside
replicates = extractReplicates(c(1:384), data)
rawData = replicates[[3]]
row.names(rawData) = rawData$sample

# Get rid of the primers that didn't work and remove their factor levels
rawData = rawData[rawData$gene.x %in% c('ATP5G3', 'B2M', 'CREB1', 'GAPDH', 'GFAP', 'GRIA2', 'HAPLN1',
                                        'LAMB1', 'OCT4', 'PROM1', 'REST'),]
rawData = droplevels(rawData)

# Extract replictes from the data set
rep1 = replicates[[1]]
rep2 = replicates[[2]]

# par(las=2)
# barchart(meanCP~gene.x,data=rawData,groups=origin.x, 
#         scales=list(x=list(rot=90,cex=0.8)), main='Raw Cp values GIC cells')
# 
# # Set some graphs
# par(las=2, mfrow=c(1,2))
# 
# # Plot the frequency of the Cp values
# hist(data$Cp, 25, main='qPCR results in raw form', col='blue', xlab='Crossing point')
# 
# # Plot correlation of the replicates
# 
# plot(rawData$Cp.x, rawData$Cp.y, main='Replicate accuracy Cp', ylab='Cp replicate 1', xlab='Cp replicate 2')
# abline(lm(rawData$Cp.x ~ rawData$Cp.y), col='red')
# summary(lm(rawData$Cp.x ~ rawData$Cp.y))
# text(50, 50, labels='R squared = 0.979')

##########################################  Generate the ddCt values by calling a custom function ######################################### 
rawData$ddCt_GAP_020_P = ddCTcalculate(geneOfInterest=rawData$gene.x, sampleOfInterest=rawData$origin.x,
                             houseKeepingGene='GAPDH', referenceSample='020_P', data=rawData)

rawData$ddCt_B2M_020_P = ddCTcalculate(geneOfInterest=rawData$gene.x, sampleOfInterest=rawData$origin.x,
                                       houseKeepingGene='B2M', referenceSample='020_P', data=rawData)

rawData$ddCt_GAP_030_P = ddCTcalculate(geneOfInterest=rawData$gene.x, sampleOfInterest=rawData$origin.x,
                                       houseKeepingGene='GAPDH', referenceSample='030_P', data=rawData)
(rawData)


#write.table(rawData, '140207_ddCt_calculations.txt', sep='\t', row.names=F)
######################################### ######################################### ######################################### 

# Subset the dataFrame by the comparisons of interest to make graphs clearer
shortLongSurvival = rawData[rawData$origin.x %in% c('020_P', '041_P'),]
primaryRecurrent = rawData[rawData$origin.x %in% c('030_P', '030a_P'),]
cd133negPos = rawData[rawData$origin.x %in% c('020_N','020_P', '030a_N', '030a_P', '041_N', '041_P'),]

# Generate ddCT values for the cd133 negative cell lines
cd133negPos$ddCt_020_N = ddCTcalculate(geneOfInterest=cd133negPos$gene.x, sampleOfInterest=cd133negPos$origin.x,
                                       houseKeepingGene='GAPDH', referenceSample='020_N', data=cd133negPos)
cd133negPos$ddCt_030_N = ddCTcalculate(geneOfInterest=cd133negPos$gene.x, sampleOfInterest=cd133negPos$origin.x,
                                       houseKeepingGene='GAPDH', referenceSample='030a_N', data=cd133negPos)
cd133negPos$ddCt_041_N = ddCTcalculate(geneOfInterest=cd133negPos$gene.x, sampleOfInterest=cd133negPos$origin.x,
                                       houseKeepingGene='GAPDH', referenceSample='041_N', data=cd133negPos)

# ######################################## Plot the ddCt values ####################################
# First compare B2M and GAPDH as house keeping genes
#################################### Try ggplot2

b2m = ggplot(data=rawData[!rawData$gene.x %in% c('GFAP'),], aes(x=origin.x, y=ddCt_B2M_020_P, fill=gene.x)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_hue(name="Gene") +      # Set legend title
    xlab("Sample") + ylab("ddCt") + # Set axis labels
    ggtitle("B2M as a house keeping gene") +  # Set title
    theme_bw(base_size=14)

gapdh = ggplot(data=rawData[!rawData$gene.x %in% c('GFAP'),], aes(x=origin.x, y=ddCt_GAP_020_P, fill=gene.x)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_hue(name="Gene") +      # Set legend title
    xlab("Sample") + ylab("ddCt") + # Set axis labels
    ggtitle("GAPDH as a house keeping gene") +  # Set title
    theme_bw(base_size=14)

multiplot(b2m, gapdh)

# Long and short term
ls1 = ggplot(data=shortLongSurvival, aes(x=origin.x, y=ddCt_B2M_020_P, fill=gene.x)) + 
           geom_bar(stat="identity", position=position_dodge(), colour="black") + 
      scale_fill_hue(name="Gene") +      # Set legend title
      xlab("Sample") + ylab("ddCt") + # Set axis labels
      ggtitle("Expression relative to long-term survivor") +  # Set title
      theme_bw(base_size=14)

ls2 = ggplot(data=shortLongSurvival[c(1:4,6:15,17:22),], aes(x=origin.x, y=ddCt_B2M_020_P, fill=gene.x)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_hue(name="Gene") +      # Set legend title
    xlab("Sample") + ylab("ddCt") + # Set axis labels
    ggtitle("Expression relative to long-term survivor") +  # Set title
    theme_bw(base_size=14)

pdf(file='140212_longShortSurv.pdf', paper='a4')
multiplot(ls1, ls2)
dev.off()

# Primary and recurrent
p1 = ggplot(data=primaryRecurrent, aes(x=origin.x, y=ddCt_GAP_030_P, fill=gene.x)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_hue(name="Gene") +      # Set legend title
    xlab("Sample") + ylab("ddCt") + # Set axis labels
    ggtitle("Expression relative to primary GIC CD133 positive") +  # Set title
    theme_bw(base_size=14)

p2 = ggplot(data=primaryRecurrent[!primaryRecurrent$gene.x %in% c('LAMB1', 'OCT4','REST'),], aes(x=origin.x, y=ddCt_GAP_030_P, fill=gene.x)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_hue(name="Gene") +      # Set legend title
    xlab("Sample") + ylab("ddCt") + # Set axis labels
    ggtitle("Expression relative to primary GIC CD133 positive") +  # Set title
    theme_bw(base_size=14)

pdf(file='140212_primaryRecurrent.pdf', paper='a4')
multiplot(p1, p2)
dev.off()

######################### Build a dataframe for individual CD133 ddCT then rbind for plotting ########
# BE CAREFUL WHEN HARD CODING ROW NAMES
# 030a
cd133_30a = rawData[c(34:55),c(1:9)]
cd133_30a$ddCt = ddCTcalculate(geneOfInterest=cd133_30a$gene.x, sampleOfInterest=cd133_30a$origin.x,
                                       houseKeepingGene='GAPDH', referenceSample='030a_N', data=cd133_30a)

#041
cd133_41 = rawData[c(56:77),c(1:9)]
cd133_41$ddCt = ddCTcalculate(geneOfInterest=cd133_41$gene.x, sampleOfInterest=cd133_41$origin.x,
                                     houseKeepingGene='GAPDH', referenceSample='041_N', data=cd133_41)

#020
cd133_20 = rawData[c(1:22),c(1:9)]
cd133_20$ddCt = ddCTcalculate(geneOfInterest=cd133_20$gene.x, sampleOfInterest=cd133_20$origin.x,
                              houseKeepingGene='GAPDH', referenceSample='020_N', data=cd133_20)

cd133negPos = rbind(cd133_20, cd133_30a, cd133_41)
# Populate a column listing whether the sample is CD133+ or CD133-
cd133negPos$cd133 = ifelse(grepl('*_N', cd133negPos$origin.x), 'negative', 'positive')
#write.table(cd133negPos, './140211_ddCtValuesCd133negPos.txt', sep='\t', row.names=T)

cd1 = ggplot(data=cd133negPos, aes(x=origin.x, y=ddCt, fill=gene.x)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_hue(name="Gene") +      # Set legend title
    xlab("Sample") + ylab("ddCt") + # Set axis labels
    ggtitle("Expression relative to CD133 negative sample") +  # Set title
    theme_bw(base_size=14)

# Take out LAMB1 to make it easier to see stuff
cd2 = ggplot(data=cd133negPos[!cd133negPos$gene.x %in% c('HAPLN1'),], aes(x=origin.x, y=ddCt, fill=gene.x)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_hue(name="Gene") +      # Set legend title
    xlab("Sample") + ylab("ddCt") + # Set axis labels
    ggtitle("Expression relative to CD133 negative sample") +  # Set title
    theme_bw(base_size=14)

pdf(file='./qcplots/ggPlot/140212_cd133NegPos.pdf', paper='a4')
multiplot(cd1, cd2)
dev.off()

# Altogerther now
multiplot(ls1, cd1, p1)
multiplot(ls2, cd2, p2)

# ################################### Test differential expression in CD133 pos neg ##############
cd133 = build_ddCTmatrix('140211_ddCtValuesCd133negPos.txt')

# Make into a numeric matrix and remove row nam
cd = as.matrix(cd133[,c(2:12)])

# Read a design matrix that I specified by hand
design = read.delim('designMatrix.txt', row.names=1)

# Need to find a group by replicate function. Two way ANOVA is the way to go
# fit = lm(cd ~ origin + cd133, design)
# x = summary(fit)
# f = lm(cd133negPos ~ origin + cd133, design)

# ################################### Basic t-test make a function one day but not 11 Feb too lazy # ################################### 
# Try a basic paired t-test
test = apply(cd, 2, pairwise.t.test, g=design$cd133, p.adjust.method='none')
y = unlist(test)
# The pvalue is contained in every 4th column
pVals = y[c(3,7,11,15,19,23,27,31,35,39,43)]
write.table(pVals, '140211_pValues.txt', sep='\t')
pVals = read.delim('140211_pValues.txt', skip=1, header=F)
colnames(pVals) = c('gene', 't_test')
pVals$p.adj = p.adjust(pVals$t_test, method='BH')

# T-test appears to give better p-values, in the end best to make it too complicated

p_valResult = rbind(cd, pVals$t_test, pVals$p.adj)
rownames(p_valResult) = c(row.names(cd), 't_test', 'p_adjust')
write.table(p_valResult, '140211_pValues.txt', sep='\t')

# ################################### gg plot make error bars for CD133 # ################################### 
# First compute summary statistics using the dataFrame generated by rbinding the individual ddCts on line 154
#cd133negPos$cd133 = ifelse(grepl('*_N', cd133negPos$origin.x), 'negative', 'positive')

# Summary statistics
cData = summariseStatistics_ddCt(cd133negPos)
cData = cbind(cData, pVals$p.adj)
colnames(cData) = c('cd133', 'gene.x','N' ,'mean', 'sd', 'se', 'adj_pVal')
cData$adj_pVal[1:11] = NA

# Add a column with stars describing if a test is significant
cData$star <- " "
cData$star[cData$adj_pVal < .05]  = "*"
cData$star[cData$adj_pVal < .01]  <- "**"
cData$star[cData$adj_pVal < .001] <- "***"
 
ggplot(cData, aes(x=gene.x, y=mean, fill=cd133)) + 
    # Define my own colors
    scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    xlab("Gene") +
    ylab("Expression normalised to GAPDH") +
    scale_fill_hue(name="CD133")+#, Legend label, use darker colors
    ggtitle("Differential expression relative to CD133 negative") +
    scale_y_continuous(breaks=0:20*4) +
    # Setting vjust to a negative number moves the asterix up a little bit to make the graph prettier
    geom_text(aes(label=star), colour="black", vjust=-2, size=10) +
    theme_bw(base_size=16)