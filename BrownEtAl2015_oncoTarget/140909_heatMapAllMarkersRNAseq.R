# Classify the RNA-seq TCGA samples provide they are different patients than Agilent
library(GSVA)
library(gplots)
library(RColorBrewer)
library(sqldf)
source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")

############################################# IO ##################################################################
setwd('~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/')

db <- dbConnect(SQLite(), dbname='~/Documents/public-datasets/cancerBrowser/deDupAgilent/coexpression.sqlite')

rnaseq = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix", row.names=1)
tcgaSigs = read.delim('~/Documents/public-datasets/TCGA/classficationSignature/131022_danFixedTCGAsignature.txt')
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

dbListTables(db)
cd133Sig = dbReadTable(db, "cd133CuttOff")
cd44Sig = dbReadTable(db, "cd44CuttOff")
cd15 = dbReadTable(db, "cd15CuttOff")
aldh1 = dbReadTable(db, "aldh1CuttOff")
itag6 = dbReadTable(db, "itag6CuttOff")
l1cam = dbReadTable(db, "l1camCuttOff")

sox2 = dbReadTable(db, "sox2CuttOff")
olig2 = dbReadTable(db, "olig2CuttOff")
ykl40 = dbReadTable(db, "ykl40CuttOff")
tubb3 = dbReadTable(db, 'tubb3CuttOff')
gfap = dbReadTable(db, 'gfapCuttOff')
id1 = dbReadTable(db, 'id1CuttOff')

pax6 = dbReadTable(db, 'pax6CuttOff')
nes = dbReadTable(db, 'nesCuttOff')
    
myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

############################################# Mung data into form for GSVA #############################################
rnaseqM = as.matrix(rnaseq)
bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), "CD15" = row.names(cd15),
               "ALDH1"=row.names(aldh1), "ITGA6"=row.names(itag6), "L1CAM"=row.names(l1cam),
               'GFAP'=row.names(gfap), 'YKL40'=row.names(ykl40), 'SOX2'=row.names(sox2), 'OLIG2'=row.names(olig2),
               'TUBB3'=row.names(tubb3), 'ID1'= row.names(id1), 'PAX6'=row.names(pax6), 'NES'=row.names(nes))

rm(cd133Sig, cd44Sig, cd15, aldh1, itag6, l1cam, olig2, gfap, ykl40, tubb3, id1, pax6, nes)

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
verhaakSubtypeAll = sort.dataframe(verhaakSubtypeAll, 'order', highFirst=F)

# The damn datatypes are not correct. Dump and read in object from file
write.table(verhaakSubtypeAll, "output.txt", sep='\t')
# write.table(verhaakSubtypeAll, "./survival/140530_liberalSignatureScores2SD.txt", sep='\t')
verhaakSubtypeAll = read.delim("output.txt", row.names=1)

subTypeHeat = as.matrix(verhaakSubtypeAll[,signatures])

# Make heat map with Veerhaak subtype
heatmap.2(t(subTypeHeat), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype and G-CIMP", #scale='row',
          Colv=as.factor(verhaakSubtypeAll$order), keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtypeAll$colours), labRow=colnames(subTypeHeat), xlab="Samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5))

############################################# Call subtypes with Agilent array daya n = 483 ##################################################################


agilentM = (as.matrix(agilentGem))
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

# Write to database
#dbWriteTable(conn = db1, name = "markerScoresAgilent", value = as.data.frame(bigResult), row.names = TRUE)
#dbDisconnect(db1)

# Merge Agilent - FACs data and clinicial data. Add Verhaak subtype
signatures = names(bigSigs)
verhaakSubtype = bindGeneExprCIMPClinical(clin, bigResult, signatures)
verhaakSubtype = sort.dataframe(verhaakSubtype, 'colours')

# The damn datatypes are not correct. Dump and read in object from file
write.table(verhaakSubtype, "output.txt", sep='\t')
# write.table(verhaakSubtype, "./survival/140603_verhaakSubtypeAgilent.txt", sep='\t')

verhaakSubtype = read.delim("output.txt", row.names=1)

subTypeHeat = as.matrix(verhaakSubtype[,signatures])

# Test ALDHA1A
verhaakSubtype$MesOther = "Mesenchymal"
verhaakSubtype$MesOther[!verhaakSubtype$GeneExp_Subtype %in% 'Mesenchymal'] = "Other"
wilcox.test(ALDH1 ~ MesOther, data=verhaakSubtype)

# Make heat map with Veerhaak subtype
heatmap.2(t(subTypeHeat), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype", 
          Colv=verhaakSubtype$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtype$colours), labRow=colnames(subTypeHeat), xlab="Aglent samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5), ylab="Marker")
