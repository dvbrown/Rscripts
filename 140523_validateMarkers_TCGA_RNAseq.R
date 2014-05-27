# Classify the RNA-seq TCGA samples provide they are different patients than Agilent
library(GSVA)
library(gplots)
library(RColorBrewer)

############################################# IO ##################################################################
setwd('~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/')
rnaseq = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix", row.names=1)
agilent = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt", row.names=1)

cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)

############################################# Mung data into form for GSVA #############################################
rnaseqM = as.matrix(rnaseq)
sigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig))

################################################# Test for enrichement #################################################

# Using ssGSEA heavily biases for CD44 subtype. Whereas for  GSVA the subtypes are more balanced
result = gsva(rnaseqM, sigs,  rnaseq=F, verbose=T, parallel.sz=1)
result = t(result$es.obs)

par(mfrow=c(2,1))
hist(result[,1], breaks='FD', main="CD133 enrichment score distribution", xlim=c(0,1))
hist(result[,2], breaks='FD', main="CD44 enrichment score distribution", xlim=c(0,1))
par(mfrow=c(1,1))

myPalette <- colorRampPalette(c("green", "black", "red"))(n = 1000)
#heatmap.2(result, cexRow=0.5, cexCol=0.9, main="ssGSEA FACS markers", 
#          scale="column", keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row")

heatmap.2(result, Colv=NA, cexRow=0.5, cexCol=0.9, main="ssGSEA FACS markers", 
          keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row")

############################################## Call the subtype #############################################
result1 = result[[2]]
subtype = ifelse(result[,1] > result[,2], "CD133", "CD44")
result1 = cbind(result, subtype)

cd133Subtype = result1[result1[,3] %in% "CD133",]
cd44Subtype = result1[result1[,3] %in% "CD44",]
length(cd133Subtype)
length(cd44Subtype)

############################################## Examine the Verhaak molecular subtypes #############################################
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

matched = intersect(row.names(clinical), colnames(rnaseq))

clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status",
                   "G_CIMP_STATUS","GeneExp_Subtype", "X_EVENT","days_to_tumor_progression", "gender")]

boundData = merge.data.frame(clin, result1)
