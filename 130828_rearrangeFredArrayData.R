#Arrange Fred's CSC microarray into a form amenable to GSEA
setwd('~/Documents/FredCSC/')
source('~/Documents//Rscripts/120704-sortDataFrame.R')
list.files(pattern='*.txt')

fdrPval = read.delim('Post-Hoc_PXR.txt',row.names=1)
avFC = read.delim('FC moyen PXR2+PXR6.txt', row.names=1)
x = read.delim('FC Pval PXR2-6.txt')

#####################################log2 transform the average fold change between PXR 2 and PXR 6###################
#can't take logarithm of a negative number
avFC$log2FC = as.numeric(abs(avFC$Fold.Change))
avFC$log2FC = log2(avFC$log2FC)
#restore the sign of the log2
avFC$log2FC = (avFC$Fold.Change / abs(avFC$Fold.Change)) * avFC$log2FC

#get rid of unecessary columns
avFC = avFC[,c(1,2,6,7,8)]
fdrPval = fdrPval[,c(1,2,8)]

#bind the fdr p-value to the log2FC. Lose all the genes for which there is no FDR available
mergeData = merge(avFC, fdrPval, by.x=c("Primary.Sequence.Name", "Sequence.Code"), 
               by.y=c("Primary.Sequence.Name", "Sequence.Code"))

#filter the data for significant measurements only. FDR of 0.1
cutoff = subset(mergeData, mergeData$Post.Hoc.P.value2.PXR2_PXR6 < 0.1)
cutoff = sort.dataframe(cutoff, 8, TRUE)

#more raw version data if wanted for GSEA
rawGSEA = sort.dataframe(avFC, 8)
rawGSEA = rawGSEA[,c(1,7,8)]

#useful columns for gsea analysis
gseaInput = cutoff[,c(1,8,9)]

setwd('./reformattedFiles/')
write.table(cutoff, '130828_filteredFDR_FredCSC.txt', sep='\t')
write.table(gseaInput, '130828_GSEAinputFilteredFDR_FredCSC.txt', sep='\t')
write.table(rawGSEA, '130828_noFDRonlyLogFC_FredCSC.txt', sep='\t')