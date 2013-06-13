#Read in survival file
setwd('~/Documents/public-datasets/TCGA/clinicalData/')
survival = read.delim('clinical_patient_gbm-sanitised.txt')
sur.sub = survival[,c(1,2,3,5,6,7,8,9,10,13,19,20,26,30)]
surv.time = sur.sub[,c(1,7)]
head(surv.time)

#strip the paitient ID to match the expression data
barcodeStrip = gsub('-F.....','-01', surv.time$bcr_followup_barcode)
numbers = surv.time$days_to_death
numbers = gsub("Not.Applicable", "NA", numbers)
numbers = gsub("Not.Available", "NA", numbers)
#numbers = gsub("[",'', numbers)

sur.strip = cbind(barcodeStrip, numbers)

sur.strip = as.data.frame(sur.strip)
colnames(sur.strip) = c('barcode', 'days_to_death')
sur.strip$days_to_death = as.numeric(sur.strip$days_to_death)
uni.sur = unique.data.frame(sur.strip)
row.names(sur.strip) = sur.strip$barcode

transpose.survival = t(sur.strip) #transpose dataframe to match the expression data
transpose.survival = as.data.frame(transpose.survival)
columns = transpose.survival[1,]
columns = as.character(columns)
colnames(transpose.survival) <- columns