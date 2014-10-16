# I aquired data from the TCGA data freeze which has combined all the data together from the various plaforms

# Methylation Patient Centric Table
# To generate DNA methylation calls for each sample per gene per overlapping platforms (HM27, HM450), we began by first collapsing multiple CpGs to one representative gene. 
# Using the associated gene expression data (organized as one gene - one expression value per sample), we merged the samples and CpG probes with gene expression data for each 
# platform. We next calculated the spearman correlation (rho) across all samples for all CpG probes for each gene to one gene expression value. For multiple CpGs for each 
# annotated gene promoter, we selected one CpG probe with the lowest correlation rho value to the associated gene expression profile to capture the most biologically 
# representative event (epigenetic silencing). This effectively reduced the number of CpG probes from N:1 to 1:1. Our data set was then reduced down to 279 samples x 9,453 
# CpG:Gene (HM27) and 74 samples x 11,040 CpG:Gene (HM450) with an overlap of 9,222 genes. Next, we assigned discrete categories based on the spearman correlation rho value 
# according to the following criteria:

# Next we assigned a ‘call’ and a confidence ‘score’ for each possible combinations (48) [3 (SNC, WNC, NNC) x 4 (CUN, CMN, VMN, IMN) x 4 (CUT, CMT, VMT, IMT)] per platform. 
# We created the following relationship for each call and score based on our interpre- tation of the most informative epigenetic event (e.g. promoter DNA hypermethylation 
# and low expression). Users should understand that the selection and criteria performed were done to the best of our knowledge at the time. 
# We felt most confident with calling epige- netically silenced events and this is reflected in the confidence score. 

# The methylation calls are as follows:
# MG: Methylation gain compared to normal ML: Methylation loss compared to normal MT: Methylated in tumor
# UT: Unmethylated in tumor
# ES: Epigenetically silenced
# UC: Unable to make call
# Methylation class confidence scores vary from 0 (no call) to 4 (high confidence).

######################################################## IO ######################################################## 
library(sqldf)
library(gplots)
library(RColorBrewer)
setwd("~/Documents/public-datasets/TCGA/methylation/")
list.files()

calls = read.delim("TCGA_GBM_dnameth_calls_20120112_ver3.txt")
scores = read.delim("TCGA_GBM_dnameth_scores_20120112_ver3.txt")
categories = read.delim("TCGA_GBM_dnameth_call_categories.txt")
data = read.delim("DNA.methylation.k6.txt")

row.names(calls) = calls$genenames.27K
calls = calls[,c(1:353)]

db <- dbConnect(SQLite(), dbname="~/Documents/public-datasets/cancerBrowser/tcgaData.sqlite")
dbListTables(db)
clinical = dbReadTable(db, "clinicalData", row.names=1)

# Retain only the ".01 cases" these arew the primary tumours.
subs = grep('*.01', row.names(clinical), value=T)
clinical = clinical[subs,]
    
# Trim the last 6 characters of the names.
row.names(clinical) = substr(row.names(clinical), 1, 15)
colnames(calls) = substr(colnames(calls), 1, 12)

# Subset the clinical data for the methylation cases that are present
clinicalCalls = clinical[colnames(calls),]
clinicalCalls = clinicalCalls[, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender", 'CDE_chemo_adjuvant_tmz', 'CDE_chemo_tmz',
                           'CDE_radiation_any', 'CDE_tmz_chemoradiation_standard', 'GeneExp_Subtype', 'G_CIMP_STATUS')]

######################################################## Assign colours to methylation calls ################################################
head(calls)
colorMatrix = apply(calls[,], c(1,2), as.character)

colorMatrix[colorMatrix == "UT"] <- 0 # Not methylated tumor
colorMatrix[colorMatrix == "MT"] <- 0.5 # Is methylated tumor
colorMatrix[colorMatrix == "ES"] <- 1 # epigenetically silenced tumor
colorMatrix[colorMatrix == "UC"] <- NA # Unable to make call
colorMatrix[colorMatrix == "ML"] <- 0.125 # methylation loss compared to normal
colorMatrix[colorMatrix == "MG"] <- 0.25 # methylation gain compared to normal

numberMatrix = apply(colorMatrix[,], c(1,2), as.numeric)
head(numberMatrix)
write.table(numberMatrix, "./141013_callsConvertNumbers.txt", sep='\t')
myPalette <- colorRampPalette(c("white", "black"))(n = 5)

######################################################## Take the top 400 genes as Brennan did ################################################

totalMeth = rowSums(numberMatrix)
cutOffMeth = quantile.default(totalMeth, probs=0.95, na.rm=T)

methylKeep = ifelse(totalMeth >= cutOffMeth, TRUE, FALSE)
methylKeep = methylKeep[methylKeep %in% TRUE]
length(methylKeep)
toBsorted = (numberMatrix[names(methylKeep),])
head(toBsorted)

# Make a heatmap where the input is a numeric representation of methylation
heatmap.2(toBsorted, cexRow=0.5,# main="Somatic mutations segrgated by marker signature", scale='none',
         # Colv=as.factor(dataSubtype$colours), key=F, trace="none", col=c('white', 'black'), 
          density.info="none", dendrogram="none", labCol='', col=myPalette,
          #ColSideColors=as.character(dataSubtype$colours), labRow=row.names(dataPresent), 
          offsetRow=c(1,1), margins=c(5,7.5)) #xlab="Samples")