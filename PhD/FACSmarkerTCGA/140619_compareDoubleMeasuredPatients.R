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

# Is there a different number of CIMP samples?
agilent = agilent[colnames(agilentGem),]
rnaSeq = rnaSeq[colnames(rnaseqGem),]
agilent$platform = "agilent"
rnaSeq$platform = 'rnaseq'
bigData = rbind(agilent, rnaSeq)
cimpTab = table(bigData$G_CIMP_STATUS, bigData$platform)
cimpTab = cimpTab[c(2,3),]
fisher.test(cimpTab)

rm(rnaseqGem, agilentGem)
################ Investigate the difference bewtween double measured patients and the rest of Agilent ###################
doublePatient = clinical[matched,]
singlePatient = clinical[!row.names(clinical) %in% matched,]

# Bind into a single dataframe
singlePatient$platform = as.factor('agilent')
doublePatient$platform = as.factor('rnaSeq')
compClinical = rbind(singlePatient, doublePatient)

# Run some tests
summary(singlePatient$days_to_death)
summary(doublePatient$days_to_death)
t.test(days_to_death ~ platform, data=compClinical)

# Measure the missingness of the data
length(singlePatient$days_to_death[is.na(singlePatient$days_to_death)])
length(doublePatient$days_to_death[is.na(doublePatient$days_to_death)])

summary(singlePatient$CDE_DxAge)
summary(doublePatient$CDE_DxAge)
t.test(CDE_DxAge ~ platform, data=compClinical) # Single patients very slightly older


t.test(doublePatient$karnofsky_performance_score, singlePatient$karnofsky_performance_score) # N.S
t.test(doublePatient$CDE_chemo_tmz_days, singlePatient$CDE_chemo_tmz_days) # N.S

# The date that the patients were treated
singleDate = as.integer(substring(singlePatient$date_of_initial_pathologic_diagnosis, 1,4))
doubleDate = as.integer(substring(doublePatient$date_of_initial_pathologic_diagnosis, 1,4))
wilcox.test(doubleDate, singleDate, conf.int=T) # double date were treated 3 years after median
median(singleDate, na.rm=T)
median(doubleDate, na.rm=T)
range(singleDate, na.rm=T)
range(doubleDate, na.rm=T)

# The theapry that the patients received. First build a contingency table
therapy = xtabs(~ CDE_therapy + platform, data=compClinical, exclude="")
# Chi squared test
summary(therapy)
therapy = table(compClinical$CDE_therapy, compClinical$platform, exclude="")
kruskal.test(list(therapy[,1], therapy[,2]))

# Construct a totals table for the fisher exact test
fisher.test(xtabs(~ CDE_therapy + platform, data=compClinical, exclude=""))
therapyDf = cbind(therapy[,1], therapy[,2])
therapyDf = cbind(therapyDf, therapyDf[,1] + therapyDf[,2])
therapyDf = rbind(therapyDf, colSums(therapyDf))
colnames(therapyDf) = c('singleMeasure', 'doubleMeasure', 'total')
fisher.test(therapyDf[1,], therapyDf[2,])

# The tissue source site (TS)
singleSource = as.integer(substring(row.names(singlePatient), 6,8))
doubleSource = as.integer(substring(row.names(doublePatient), 6,8))
wilcox.test(singleSource, doubleSource, conf.int=F)

# The molecular subtype is not significant
molSub = xtabs(~ GeneExp_Subtype + platform, data= compClinical, exclude="")
summary(molSub)

# G-CIMP is not significant
fisher.test(xtabs(~ G_CIMP_STATUS + platform, data= compClinical, exclude=""))

# Set up multiple fisher tests to investigate clinical parameters
fisher.test(xtabs(~ gender + platform, data=compClinical))# N.S
fisher.test(xtabs(~ CDE_alk_chemoradiation_standard + platform, data=compClinical))# N.S

fisher.test(xtabs(~ CDE_chemo_tmz + platform, data=compClinical, exclude = "")) # 0.011

fisher.test(xtabs(~ CDE_chemo_tmz_long + platform, data=compClinical)) # N.S
fisher.test(xtabs(~ CDE_radiation_standard + platform, data=compClinical)) # N.S
fisher.test(xtabs(~ CDE_radiation_adjuvant + platform, data=compClinical)) # N.S
fisher.test(xtabs(~ CDE_radiation_adjuvant_standard + platform, data=compClinical)) # NS
t.test(CDE_chemo_alk_days ~ platform, data=compClinical) # N.S

fisher.test(xtabs(~ pretreatment_history + platform, data=compClinical, exclude = "")) # NS    
fisher.test(xtabs(~ CDE_chemo_adjuvant_tmz + platform, data=compClinical, exclude = "")) # Highly significant    
fisher.test(xtabs(~ CDE_chemo_adjuvant_alk + platform, data=compClinical, exclude = "")) # almost significant        
fisher.test(xtabs(~ additional_surgery_locoregional_procedure + platform, data=compClinical, exclude = "")) # barely significant

# Read in the hand made summary table and compute FDR
# clinTests = read.delim('~/Documents/public-datasets/cancerBrowser/compareDataSourcePlots/140702_statsTestingClinical_transposed.txt')
# clinTests$fdr = p.adjust(clinTests$p.value)#, method='fdr')
# write.table(clinTests, '~/Documents/public-datasets/cancerBrowser/compareDataSourcePlots/140702_statsTestingClinical_fdr.txt', sep='\t', row.names=F)

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
plot(surFitAgilent, main='TCGA GBM cohort classified by FACS marker signature Agilent array', ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),
     xlim=c(0,1600), 
     cex=1.75, conf.int=F, lwd=1.33)
     
legend('topright', c('CD133', 'CD44'),
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

plot(surFitRNA, main='TCGA GBM cohort classified by FACS marker signature RNAseq', ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),
     xlim=c(0,1600), 
     cex=1.75, conf.int=F, lwd=1.33)

legend('topright', c('CD133', 'CD44'),
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

#test for a difference between curves
testA = surv_test(survAgilent~boundDataA$subtype)#, subset=!boundData$subtype %in% "intermediate")
testA
testR = surv_test(survRNA~boundDataR$subtype)#, subset=!boundData$subtype %in% "intermediate")
testR