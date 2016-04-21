# Classify the RNA-seq TCGA samples provide they are different patients than Agilent
library(GSVA)
library(gplots)
library(RColorBrewer)
source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")

############################################# IO ##################################################################
setwd('~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/')
agilent = read.delim('~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt', row.names=1)# agilent = read.delim('~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt', row.names=1)
tcgaSigs = read.delim('~/Documents/public-datasets/TCGA/classficationSignature/131022_danFixedTCGAsignature.txt')

clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)
cd15 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/CD15/140528_cd15Cutoff.txt", row.names=1)
aldh1 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ALDH1/140528_ALDH1Cutoff.txt", row.names=1)
itag6 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ITGA6//140528_ITGA6Cutoff.txt", row.names=1)
l1cam = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/L1CAM/140528_L1CAMCutoff.txt", row.names=1)

myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

############################################# Call subtypes with Agilent array daya n = 483 ##################################################################
agilentM = t(as.matrix(agilent))
# Change NAs to 0 as GSVA doesn't like it
agilentM[is.na(agilentM)] <- 0

bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), "CD15" = row.names(cd15),
               "ALDH1"=row.names(aldh1), "ITGA6"=row.names(itag6), "L1CAM"=row.names(l1cam))
rm(cd133Sig, cd44Sig, cd15, aldh1, itag6, l1cam)

# Extract the clinical data for the Agilent patients
matched = intersect(row.names(clinical), colnames(agilentM))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status",
                           "G_CIMP_STATUS","GeneExp_Subtype", "X_EVENT","days_to_tumor_progression", "gender")]

# Call the subtypes with GSVA
bigResult = gsva(agilentM, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1)
bigResult = t(bigResult$es.obs)

# Merge Agilent - FACs data and clinicial data. Add Verhaak subtype
signatures = names(bigSigs)
verhaakSubtype = bindGeneExprCIMPClinical(clin, bigResult, signatures)
verhaakSubtype = sort.dataframe(verhaakSubtype, 'GeneExp_Subtype')

# The damn datatypes are not correct. Dump and read in object from file
write.table(verhaakSubtype, "output.txt", sep='\t')
# write.table(verhaakSubtype, "./survival/140603_verhaakSubtypeAgilent.txt", sep='\t')

verhaakSubtype = read.delim("output.txt", row.names=1)
verhaakSubtype = verhaakSubtype[!verhaakSubtype$colours %in% 'black',]
subTypeHeat = as.matrix(verhaakSubtype[,signatures])

# Make heat map with Veerhaak subtype
heatmap.2(t(subTypeHeat), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype", 
          Colv=verhaakSubtype$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtype$colours), labRow=colnames(subTypeHeat), xlab="Agilent samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5))

# Make heat map with only 3 markers
heatmap.2(t(subTypeHeat[,c(1:3)]), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype and G-CIMP", #scale='row',
          Colv=verhaakSubtype$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtype$colours), labRow=colnames(subTypeHeat), xlab="Agilent samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5))

# Make heat map with only 3 mRNAs
markers = as.matrix(agilentM[c('CD44', 'FUT4', 'PROM1'),])
markers = markers[,row.names(verhaakSubtype)]

heatmap.2(markers, cexRow=1.5, main="Enrichment of FACS marker mRNAs\n in Molecular Subtype and G-CIMP", scale='row',
          Colv=verhaakSubtype$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtype$colours), labRow=colnames(subTypeHeat), xlab="Agilent samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5))

############################################# Us the non default enrichment statistic ##################################################################
# Call the subtypes with GSVA
bigResult2 = gsva(agilentM, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1, mx.diff=FALSE)
bigResult2 = (bigResult2$es.obs)

# Merge Agilent - FACs data and clinicial data. Add Verhaak subtype
signatures = names(bigSigs)
verhaakSubtype = bindGeneExprCIMPClinical(clin, bigResult, signatures)
verhaakSubtype = sort.dataframe(verhaakSubtype, 'GeneExp_Subtype')

# The damn datatypes are not correct. Dump and read in object from file
write.table(verhaakSubtype, "output.txt", sep='\t')
# write.table(verhaakSubtype, "./survival/140603_verhaakSubtypeAgilent.txt", sep='\t')

verhaakSubtype = read.delim("output.txt", row.names=1)
verhaakSubtype = verhaakSubtype[!verhaakSubtype$colours %in% 'black',]
subTypeHeat = as.matrix(verhaakSubtype[,signatures])

# Make an object to plot density gram
lattPlot = data.frame((cbind(bigResult[,"CD133"], bigResult2["CD133",])))#, 

par(mfrow=c(2,1))
hist(bigResult[,"CD133"])
hist(bigResult2["CD133",])