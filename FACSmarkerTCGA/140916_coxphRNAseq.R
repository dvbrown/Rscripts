library(survival)
library(coin)
library(ggplot2)

source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")

############################################## IO and general munging #############################################
# Load the signature
verhaakSignature = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140530_liberalSignatureScores2SD.txt", row.names=1)

verhaakSignature = verhaakSignature[,c("CD133","CD44","GeneExp_Subtype","G_CIMP_STATUS", 'colours')]
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

# Call the FACS subtype
verhaakSubtypeCall = callMarkerSubtype(verhaakSignature, 0, 0)

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), row.names(verhaakSubtypeCall))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender", 'CDE_chemo_adjuvant_tmz', 'CDE_chemo_tmz',
                           'CDE_radiation_any', 'CDE_tmz_chemoradiation_standard', 'GeneExp_Subtype')]

############################################## bind the clinical and subtyping info together #############################################

boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)

# Fix up the GCIMP to only have true and false
boundData$G_CIMP_STATUS = as.character(boundData$G_CIMP_STATUS)
boundData$G_CIMP_STATUS[boundData$G_CIMP_STATUS %in% 'G-CIMP'] = TRUE
boundData$G_CIMP_STATUS[!boundData$G_CIMP_STATUS %in% 'TRUE'] = FALSE
boundData$G_CIMP_STATUS = as.factor(boundData$G_CIMP_STATUS)

############################################# Analysing the data for survival ##################################
data.surv = Surv(boundData$CDE_survival_time, event=boundData$X_EVENT)
coxPH = coxph(data.surv ~  subtype +  CDE_DxAge + CDE_chemo_tmz+  CDE_radiation_any + gender, 
              data=boundData, na.action="na.omit")
summary(coxPH)