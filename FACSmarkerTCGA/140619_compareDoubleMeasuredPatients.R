# Comapre the same patients that were RNAseq'd and Agilented
library(GSVA)
library(gplots)
library(RColorBrewer)
library(survival)
library(coin)
source("~/Documents/Rscripts/120704-sortDataFrame.R")
source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")
setwd('~/Documents/public-datasets/cancerBrowser/dedupAgilent/')

myPalette <- colorRampPalette(c("green", "black", "red"))(n = 1000)

clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)
# subset the clinical data for only primary GBMs
clinical = clinical[clinical$histological_type %in% c('Untreated primary (de novo) GBM', 'Glioblastoma Multiforme (GBM)'),]

rnaSeq = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)
rnaseqGem = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix", row.names=1)

agilent = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_G4502A_07_2-2014-05-02/clinical_dataDots.txt", row.names=1)
agilentGem = read.delim('~/Documents/public-datasets/cancerBrowser/TCGA_GBM_G4502A_07_2-2014-05-02/genomicMatrix', row.names=1)

matched = intersect(colnames(rnaseqGem), colnames(agilentGem))

################ Investigate the difference bewtween double measured patients and the rest of Agilent ###################
doublePatient = clinical[matched,]
singlePatient = clinical[!row.names(clinical) %in% matched,]

# Bind into a single dataframe
singlePatient$platform = as.factor('agilent')
doublePatient$platform = as.factor('rnaSeq')
compClinical = rbind(singlePatient, doublePatient)

# Run some tests
t.test(doublePatient$days_to_death, singlePatient$days_to_death) # Sig large
t.test(doublePatient$CDE_DxAge, singlePatient$CDE_DxAge) # Single patients very slightly older
t.test(doublePatient$karnofsky_performance_score, singlePatient$karnofsky_performance_score) # N.S
t.test(doublePatient$CDE_chemo_tmz_days, singlePatient$CDE_chemo_tmz_days) # N.S

# The date that the patients were treated
singleDate = as.integer(substring(singlePatient$date_of_initial_pathologic_diagnosis, 1,4))
doubleDate = as.integer(substring(doublePatient$date_of_initial_pathologic_diagnosis, 1,4))
wilcox.test(doubleDate, singleDate, conf.int=T) # double date were treated 3 years after median
median(singleDate, na.rm=T)
median(doubleDate, na.rm=T)

# The theapry that the patients received. First build a contingency table
therapy = xtabs(~ CDE_therapy + platform, data=compClinical, exclude="")
# Chi squared test
summary(therapy)
therapy = table(compClinical$CDE_therapy, compClinical$platform, exclude="")
kruskal.test(list(therapy[,1], therapy[,2]))
prop.table(therapy[,1], therapy[,2])

# The tissue source site (TS)
singleSource = as.integer(substring(row.names(singlePatient), 6,8))
doubleSource = as.integer(substring(row.names(doublePatient), 6,8))
wilcox.test(singleSource, doubleSource, conf.int=F)

################ Mung the GEMs into shape for GSVA #############################

rnaseqGemM = rnaseqGem[,matched]
agilentGemM = agilentGem[,matched]

par(mfrow=c(2,1))
boxplot(agilentGemM, main='Matched Agilent')
boxplot(rnaseqGemM, main='Matched RNAseq')

#### The signature computation ####

cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)

agilentM = as.matrix(agilentGemM)
rnaSeqM = as.matrix(rnaseqGemM)

# Change NAs to 0 as GSVA doesn't like it
agilentM[is.na(agilentM)] <- 0

#Fix the row.names of the clinical names to atch dfata
#row.names(agilentM) = paste(row.names(agilentM), '.01', sep='')

bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig))

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), colnames(agilentM))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status",
                           "G_CIMP_STATUS","GeneExp_Subtype", "X_EVENT","days_to_tumor_progression", "gender")]

# Call the subtypes with GSVA
resultAgilent = gsva(agilentM, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1)
resultAgilent = t(resultAgilent$es.obs)
resultRNAseq = gsva(rnaSeqM, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1)
resultRNAseq = t(resultRNAseq$es.obs)

# Merge Agilent - FACs data and clinicial data. Add Verhaak subtype
signatures = names(bigSigs)

subtypeAgilent = bindGeneExprCIMPClinical(clin, resultAgilent, signatures)
subtypeAgilent = sort.dataframe(subtypeAgilent, 'colours')
write.table(subtypeAgilent, "output.txt", sep='\t')
subtypeAgilent = read.delim("output.txt", row.names=1)
subTypeHeatAgilent = as.matrix(subtypeAgilent[,signatures])

subtypeRNA = bindGeneExprCIMPClinical(clin, resultRNAseq, signatures)
subtypeRNA = sort.dataframe(subtypeRNA, 'colours')
write.table(subtypeRNA, "output.txt", sep='\t')
subtypeRNA = read.delim("output.txt", row.names=1)
subTypeHeatRNA = as.matrix(subtypeRNA[,signatures])

# Make heat map with Veerhaak subtype
par(mfrow=c(2,1))

heatmap.2(t(subTypeHeatAgilent), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype", 
          Colv=subtypeAgilent$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(subtypeAgilent$colours), labRow=colnames(subtypeAgilent), xlab="Aglent samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5), ylab="Marker")

heatmap.2(t(subTypeHeatRNA), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype", 
          Colv=subtypeRNA$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(subtypeRNA$colours), labRow=colnames(subtypeRNA), xlab="RNA samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5), ylab="Marker")

subtypeRNACall = callMarkerSubtype(subtypeRNA, 0, 0)
subtypeAgilentCall = callMarkerSubtype(subtypeAgilent, 0, 0)

############################################## Survival analysis #############################################
boundDataA = merge.data.frame(clin, subtypeAgilentCall, by.x="row.names", by.y="row.names")
boundDataA = sort.dataframe(boundDataA, "subtype")
row.names(boundDataA) = boundDataA$Row.names
boundDataA$subtype = as.factor(boundDataA$subtype)
boundDataA$gender = as.factor(boundDataA$gender)

boundDataR = merge.data.frame(clin, subtypeRNACall, by.x="row.names", by.y="row.names")
boundDataR = sort.dataframe(boundDataR, "subtype")
row.names(boundDataR) = boundDataR$Row.names
boundDataR$subtype = as.factor(boundDataR$subtype)
boundDataR$gender = as.factor(boundDataR$gender)

############################################# Analysing the data for survival ##################################

#generate the survival object and plot a Kaplan-Meier
survAgilent = Surv(boundDataA$CDE_survival_time, event=boundDataA$X_EVENT)
survRNA = Surv(boundDataR$CDE_survival_time, event=boundDataR$X_EVENT)

surFitAgilent = survfit(survAgilent~subtype, boundDataA)
surFitRNA = survfit(survRNA~subtype, boundDataR)

par(mfrow=c(2,1))
plot(surFitAgilent, main='FACS marker coexpression signature in \nGlioblastoma multiforme Agilent duplicates',
     ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),#'green'),
     xlim=c(0,750), cex=1.75, conf.int=F, lwd=1.5)
legend('topright', c('CD133', 'CD44'),# 'Intermediate'), 
       col=c("red",'blue'),#'green'),
       lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

plot(surFitRNA, main='FACS marker coexpression signature in \nGlioblastoma multiforme RNAseq duplicates',
     ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),#'green'),
     xlim=c(0,750), cex=1.75, conf.int=F, lwd=1.5)
legend('topright', c('CD133', 'CD44'),# 'Intermediate'), 
       col=c("red",'blue'),#'green'),
       lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

#test for a difference between curves
testA = surv_test(survAgilent~boundDataA$subtype)#, subset=!boundData$subtype %in% "intermediate")
testA
testR = surv_test(survRNA~boundDataR$subtype)#, subset=!boundData$subtype %in% "intermediate")
testR