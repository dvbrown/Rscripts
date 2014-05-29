# Classify the RNA-seq TCGA samples provide they are different patients than Agilent
library(GSVA)
library(gplots)
library(RColorBrewer)
source("~/Documents/Rscripts/120704-sortDataFrame.R")

bindGeneExprClinical <- function (clinicalData, subtypedGeneExpression, signatures) {
    # Merges clinical and FACS marker subtyped gene expression information and annotate a color based on Verhaak subtype
    # signatures is a character vector of the signature names
    boundData = merge.data.frame(clinicalData, subtypedGeneExpression, by.x="row.names", by.y="row.names")
    verhaakSubtype = boundData[,c(signatures, "GeneExp_Subtype")]
    verhaakSubtype$colours = "black"
    verhaakSubtype$colours[verhaakSubtype$GeneExp_Subtype == "Proneural"] = "red"
    verhaakSubtype$colours[verhaakSubtype$GeneExp_Subtype == "Neural"] = "green"
    verhaakSubtype$colours[verhaakSubtype$GeneExp_Subtype == "Classical"] = "blue"
    verhaakSubtype$colours[boundData$GeneExp_Subtype == "Mesenchymal"] = "orange"
    return (verhaakSubtype)
}

bindGeneExprCIMPClinical <- function (clinicalData, subtypedGeneExpression, signatures) {
    # Merges clinical and FACS marker subtyped gene expression information and annotate a color based on Verhaak subtype
    # signatures is a character vector of the signature names
    boundData = merge.data.frame(clinicalData, subtypedGeneExpression, by.x="row.names", by.y="row.names")
    verhaakSubtype = boundData[,c(signatures, "GeneExp_Subtype", "G_CIMP_STATUS")]
    verhaakSubtype$colours = "black"
    verhaakSubtype$colours[verhaakSubtype$GeneExp_Subtype == "Proneural"] = "red"
    verhaakSubtype$colours[verhaakSubtype$GeneExp_Subtype == "Neural"] = "green"
    verhaakSubtype$colours[verhaakSubtype$GeneExp_Subtype == "Classical"] = "blue"
    verhaakSubtype$colours[verhaakSubtype$GeneExp_Subtype == "Mesenchymal"] = "orange"
    verhaakSubtype$colours[verhaakSubtype$G_CIMP_STATUS == "G-CIMP"] = "pink"
    return (verhaakSubtype)
}

############################################# IO ##################################################################
setwd('~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/')
rnaseq = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix", row.names=1)
# agilent = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt", row.names=1)
tcgaSigs = read.delim('~/Documents/public-datasets/TCGA/classficationSignature/131022_danFixedTCGAsignature.txt')

clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)
cd15 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/CD15/140528_cd15Cutoff.txt", row.names=1)
aldh1 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ALDH1/140528_ALDH1Cutoff.txt", row.names=1)
itag6 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ITGA6//140528_ITGA6Cutoff.txt", row.names=1)
l1cam = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/L1CAM/140528_L1CAMCutoff.txt", row.names=1)

myPalette <- colorRampPalette(c("green", "black", "red"))(n = 1000)

############################################# Mung data into form for GSVA #############################################
rnaseqM = as.matrix(rnaseq)
sigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig))

################################################# Test for enrichment #################################################

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

#heatmap.2(result, cexRow=0.5, cexCol=0.9, main="ssGSEA FACS markers", 
#          scale="column", keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row")

heatmap.2(result, Colv=NA, cexRow=0.5, cexCol=0.9, main="ssGSEA FACS markers", 
          keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row")

############################################## Examine the Verhaak molecular subtypes for CD133 and CD44 #############################################
# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), colnames(rnaseq))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status",
                   "G_CIMP_STATUS","GeneExp_Subtype", "X_EVENT","days_to_tumor_progression", "gender")]

# Merge RNAseq - FACs data and clinicial data. Add Verhaak subtype
verhaakSubtype = bindGeneExprClinical(clin, result1, c("CD133", "CD44"))
verhaakSubtype = sort.dataframe(verhaakSubtype, 4)

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

############################################## Heatmap all signatures with Verhaak molecular subtypes #############################################
bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), "CD15" = row.names(cd15),
               "ALDH1"=row.names(aldh1), "ITGA6"=row.names(itag6), "L1CAM"=row.names(l1cam))

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
verhaakSubtypeAll = read.delim("output.txt", row.names=1)

subTypeHeat = as.matrix(verhaakSubtypeAll[,signatures])

# Make heat map with Veerhaak subtype
heatmap.2(t(subTypeHeat), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype and G-CIMP", 
          Colv=verhaakSubtypeAll$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtypeAll$colours), labRow=colnames(subTypeHeat), xlab="Samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5))

############################################## Test for differential CD133 expression of proneural and G-CIMP #############################################

subtype = ifelse(verhaakSubtypeAll[,"CD133"] > verhaakSubtypeAll[,"CD44"], "CD133", "CD44")
verhaakFACSSubtypeAll = cbind(verhaakSubtypeAll, subtype)

# Subset only the proneural and G-CIMP cases
verhaakFACSSubtypeAll = verhaakFACSSubtypeAll[verhaakFACSSubtypeAll$GeneExp_Subtyp %in% "Proneural" | verhaakFACSSubtypeAll$G_CIMP_STATUS %in% "G-CIMP",]
verhaakFACSSubtypeAll$GeneExp_Subtype[verhaakFACSSubtypeAll$G_CIMP_STATUS == "G-CIMP"] = "G-CIMP"
verhaakFACSSubtypeAll = droplevels(verhaakFACSSubtypeAll)

# t.test for differential expression only G-CIMP and CD133 cases
t.test(CD133 ~ GeneExp_Subtype, verhaakFACSSubtypeAll)