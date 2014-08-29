# Classify the RNA-seq TCGA samples provide they are different patients than Agilent
library(GSVA)
library(gplots)
library(RColorBrewer)
source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")

############################################# IO ##################################################################
setwd('~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/')
rnaseq = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix", row.names=1)
# agilent = read.delim('~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt', row.names=1)
tcgaSigs = read.delim('~/Documents/public-datasets/TCGA/classficationSignature/131022_danFixedTCGAsignature.txt')

clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)
cd15 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/CD15/140528_cd15Cutoff.txt", row.names=1)
aldh1 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ALDH1/140528_ALDH1Cutoff.txt", row.names=1)
itag6 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ITGA6//140528_ITGA6Cutoff.txt", row.names=1)
l1cam = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/L1CAM/140528_L1CAMCutoff.txt", row.names=1)

myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

############################################# Mung data into form for GSVA #############################################
rnaseqM = as.matrix(rnaseq)
bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), "CD15" = row.names(cd15),
               "ALDH1"=row.names(aldh1), "ITGA6"=row.names(itag6), "L1CAM"=row.names(l1cam))
rm(cd133Sig, cd44Sig, cd15, aldh1, itag6, l1cam)

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), colnames(rnaseq))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status",
                           "G_CIMP_STATUS","GeneExp_Subtype", "X_EVENT","days_to_tumor_progression", "gender")]

############################################## Heatmap all signatures with Verhaak molecular subtypes #############################################

# Using ssGSEA heavily biases for CD44 subtype. Whereas for  GSVA the subtypes are more balanced
bigResult = gsva(rnaseqM, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1)
bigResult = t(bigResult$es.obs)

# Merge RNAseq - FACs data and clinicial data. Add Verhaak subtype
signatures = names(bigSigs)
verhaakSubtype = bindGeneExprClinical(clin, bigResult, signatures)
verhaakSubtype = sort.dataframe(verhaakSubtype, 'colours')

# The damn datatypes are not correct. Dump and read in object from file
write.table(verhaakSubtype, "output.txt", sep='\t')
verhaakSubtype = read.delim("output.txt", row.names=1)

subTypeHeat = as.matrix(verhaakSubtype[,signatures])

# Make heat map with Veerhaak subtype
heatmap.2(t(subTypeHeat), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype", 
          Colv=verhaakSubtype$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtype$colours), labRow=colnames(subTypeHeat), xlab="Samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5))

####################################### Heatmap all signatures with Verhaak molecular subtypes and G-CIMP #############################################

verhaakSubtypeAll = bindGeneExprCIMPClinical(clin, bigResult, signatures)
verhaakSubtypeAll = sort.dataframe(verhaakSubtypeAll, 'colours')

# The damn datatypes are not correct. Dump and read in object from file
write.table(verhaakSubtypeAll, "output.txt", sep='\t')
# write.table(verhaakSubtypeAll, "./survival/140530_liberalSignatureScores2SD.txt", sep='\t')
verhaakSubtypeAll = read.delim("output.txt", row.names=1)

subTypeHeat = as.matrix(verhaakSubtypeAll[,signatures])

# Make heat map with Veerhaak subtype
heatmap.2(t(subTypeHeat), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype and G-CIMP", #scale='row',
          Colv=verhaakSubtypeAll$G_CIMP_STATUS, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtypeAll$colours), labRow=colnames(subTypeHeat), xlab="Samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5))

# Make heat map with only 3 markers
heatmap.2(t(subTypeHeat[,c(1:3)]), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype and G-CIMP", #scale='row',
          Colv=verhaakSubtypeAll$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtypeAll$colours), labRow=colnames(subTypeHeat), xlab="Samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5))

# Make heat map with only 3 mRNAs
markers = as.matrix(rnaseq[c('CD44', 'FUT4', 'PROM1'),])

heatmap.2(markers, cexRow=1.5, main="Enrichment of FACS marker mRNAs\n in Molecular Subtype and G-CIMP", #scale='row',
          Colv=verhaakSubtypeAll$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtypeAll$colours), labRow=row.names(markers), xlab="Samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5))


################################################# Test for enrichment of only CD133 and CD44 using the 3SD signature #################################################

cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140529_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140529_cd44Cutoff.txt", row.names=1)

sigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig))

# Using ssGSEA heavily biases for CD44 subtype. Whereas for  GSVA the subtypes are more balanced
result = gsva(rnaseqM, sigs,  rnaseq=F, verbose=T, parallel.sz=1)
result = t(result$es.obs)

# Compare ssGSEA
# result = gsva(rnaseqM, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1, method="ssgsea")
# result = t(result)

par(mfrow=c(2,1))
hist(result[,1], breaks='FD', main="CD133 enrichment score distribution", xlim=c(0,1))
hist(result[,2], breaks='FD', main="CD44 enrichment score distribution", xlim=c(0,1))
par(mfrow=c(1,1))

heatmap.2(result, Colv=NA, cexRow=0.5, cexCol=0.9, main="ssGSEA FACS markers", 
          keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row")

############################################## Make a heatmap for the 3SD signature #############################################

# Merge RNAseq - FACs data and clinicial data. Add Verhaak subtype
verhaakSubtype = bindGeneExprClinical(clin, result, c("CD133", "CD44"))
verhaakSubtype = sort.dataframe(verhaakSubtype, 'colours')

# The damn datatypes are not correct. Dump and read in object from file
write.table(verhaakSubtype, "output.txt", sep='\t')
verhaakSubtype = read.delim("output.txt", row.names=1)

subTypeHeat = as.matrix(verhaakSubtype[,c("CD133","CD44")])

# Make heat map with Veerhaak subtype
heatmap.2(t(subTypeHeat), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype", Colv=verhaakSubtype$colours,
          keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="none", ColSideColors=as.character(verhaakSubtype$colours),
          labRow=c("CD133", "CD44"), xlab="Samples", labCol=NA, offsetRow=c(1,1), margins=c(2,7))

############################################## Test for mutual exclusivity of proneural and mesenchymal #############################################
subtype = ifelse(verhaakSubtype[,"CD133"] > verhaakSubtype[,"CD44"], "CD133", "CD44")
verhaakFACSSubtype = cbind(verhaakSubtype, subtype)

# Subset only the proneural and Mesenchymal cases
verhaakFACSSubtype = verhaakFACSSubtype[verhaakFACSSubtype$GeneExp_Subtyp %in% "Proneural" | verhaakFACSSubtype$GeneExp_Subtyp %in% "Mesenchymal",]
verhaakFACSSubtype = droplevels(verhaakFACSSubtype)

# Build a contingency table
contingency = table(verhaakFACSSubtype[,c("GeneExp_Subtype", "subtype")])
fisher.test(contingency)

############################################## Test for differential CD133 expression of proneural and G-CIMP #############################################
# Add the G-CIMP annotation
verhaakSubtype = bindGeneExprCIMPClinical(clin, result, c("CD133", "CD44"))
verhaakSubtype = sort.dataframe(verhaakSubtype, 'G_CIMP_STATUS')
# write.table(verhaakSubtype, "./survival/140529_verhaakSubtypeCD133_scores", sep='\t')

# # The damn datatypes are not correct. Dump and read in object from file
# write.table(verhaakSubtype, "output.txt", sep='\t')
# verhaakSubtype = read.delim("output.txt", row.names=1)

# Subset only the proneural and G-CIMP cases
verhaakSubtype = verhaakSubtype[verhaakSubtype$GeneExp_Subtyp %in% "Proneural" | verhaakSubtype$G_CIMP_STATUS %in% "G-CIMP",]
verhaakSubtype$GeneExp_Subtype[verhaakSubtype$G_CIMP_STATUS == "G-CIMP"] = "G-CIMP"
verhaakSubtype = droplevels(verhaakSubtype)

# t.test for differential expression only G-CIMP and CD133 cases for the 3SD cutoff
t.test(CD133 ~ GeneExp_Subtype, verhaakSubtype)

# Make heat map with Veerhaak subtype
heatmap.2(t(subTypeHeat), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype and G-CIMP", 
          Colv=verhaakSubtypeAll$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtypeAll$colours), labRow=colnames(subTypeHeat), xlab="Samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5))



############################################# Call subtypes with Agilent array daya n = 483 ##################################################################
agilentM = t(as.matrix(agilent))
# Change NAs to 0 as GSVA doesn't like it
agilentM[is.na(agilentM)] <- 0

bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), "CD15" = row.names(cd15),
               "ALDH1"=row.names(aldh1), "ITGA6"=row.names(itag6), "L1CAM"=row.names(l1cam))

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
verhaakSubtype = sort.dataframe(verhaakSubtype, 'colours')

# The damn datatypes are not correct. Dump and read in object from file
write.table(verhaakSubtype, "output.txt", sep='\t')
# write.table(verhaakSubtype, "./survival/140603_verhaakSubtypeAgilent.txt", sep='\t')

verhaakSubtype = read.delim("output.txt", row.names=1)

subTypeHeat = as.matrix(verhaakSubtype[,signatures])

# Make heat map with Veerhaak subtype
heatmap.2(t(subTypeHeat), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype", 
          Colv=verhaakSubtype$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtype$colours), labRow=colnames(subTypeHeat), xlab="Aglent samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5), ylab="Marker")

############################################# Try the 3SD signature gsva ##################################
# Classify the RNA-seq TCGA samples provide they are different patients than Agilent
library(GSVA)
library(gplots)
library(RColorBrewer)
source("~/Documents/Rscripts/120704-sortDataFrame.R")
source("~/Documents/Rscripts/140508_coexpressionFunctions.R")

tcgaSigs = read.delim('~/Documents/public-datasets/TCGA/classficationSignature/131022_danFixedTCGAsignature.txt')

clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140529_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140529_cd44Cutoff.txt", row.names=1)
agilent = read.delim('~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt', row.names=1)
agilentM = t(as.matrix(agilent))
# Change NAs to 0 as GSVA doesn't like it
agilentM[is.na(agilentM)] <- 0

sigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig))

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), colnames(agilentM))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status",
                           "G_CIMP_STATUS","GeneExp_Subtype", "X_EVENT","days_to_tumor_progression", "gender")]

# Call the subtypes with GSVA
bigResult = gsva(agilentM, sigs,  rnaseq=F, verbose=T, parallel.sz=1, method="ssgsea")
bigResult = t(bigResult$es.obs)

# Merge Agilent - FACs data and clinicial data. Add Verhaak subtype
signatures = names(sigs)
verhaakSubtype = bindGeneExprCIMPClinical(clin, bigResult, signatures)
verhaakSubtype = sort.dataframe(verhaakSubtype, 'colours')

# The damn datatypes are not correct. Dump and read in object from file
write.table(verhaakSubtype, "output.txt", sep='\t')
# write.table(verhaakSubtype, "./survival/140606_verhaakSubtypeAgilent_3sd.txt", sep='\t')

verhaakSubtype = read.delim("output.txt", row.names=1)

subTypeHeat = as.matrix(verhaakSubtype[,signatures])

# Make heat map with Veerhaak subtype
heatmap.2(t(subTypeHeat), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype", 
          Colv=verhaakSubtype$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtype$colours), labRow=colnames(subTypeHeat), xlab="Aglent samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5), ylab="Marker")
