#Take the CREB signature and do a Cox prop regression on Gravandeel et al
library(survival)
setwd('~/Documents/CREB/gliomaMicroarray/Gravendeel/analysis/')
source('~/Documents/Rscripts/121112_subSetTCGAfunctions.R')

extractDataClinical = function(clinicalData) { #A function to get useful clinical data and match rownames to expression
  sur = clinicalData[,,]
  barcode = as.character(clinicalData[,1])
  barcode = gsub('$', '.CEL', barcode) #append .CEL to the first row to match the array filenames
  sur = cbind(barcode, sur)
  row.names(sur) = sur$barcode
  return (sur)
}

bindGEMsurvival = function(geneExpressionMatrix, clinical) {
  #bind the computed gene signature score to clinical data for analysis
  set = intersect(row.names(geneExpressionMatrix), row.names(clinical)) #get the overlap
  data = cbind(geneExpressionMatrix[set,], clinical[set,]) #bind the 2 datasets together based on overlap
  colnames(data) = c(sig, 'barcode', 'filename','number', 'sex', 'grade','status','survivalMonth', 'survivalYear','grade.code','age')
  data$survivalMonth = as.numeric(data$survivalMonth) 
  data$survivalYear = as.numeric(data$survivalYear)
  data$age = as.numeric(data$age)
  data$grade.code=as.factor(data$grade.code)
  return(data)
}

#prepare the data for multiple regression
data = read.delim('130123_zScoreGravendeelTransformeded.txt', sep='\t',row.names=1)
#CREB signature change this!
signature = read.delim('~/Documents/CREB/CREBsignatures/LembergerRamsayOverlap/130123_LembergerCREBkoOnlyTop50.txt')
signature = as.character(signature[,1]) #the score column ranks the number of sites in the given promoter

clinical = read.delim('~/Documents/CREB/gliomaMicroarray/Gravendeel/analysis/designV2.txt')
clinical = extractDataClinical(clinical)
sig = intersect(row.names(data),signature)
subsetGrav = data[sig,]
data.t = t(subsetGrav)
finalData = bindGEMsurvival(data.t, clinical)

#got to comment out the geneList for reuse
#LembergerOnlyTop50 = AVPI1+C1QA+C1QB+C3AR1+CCDC116+CD14+COX17+CRLS1+CST7+CTSC+CTSE+CXCL13+DMRT1+DNM2+EIF6+FCER1G+GSTA3+ITGAX+LPHN2+LRPAP1+MCM5+MORN1+NKRF+NRG4+OSMR+PAFAH1B2+PCYOX1+PRAME+RAB11A+RAP2A+RICTOR+SRSF7+TGFBR2+THTPA+TLR2+TMEM117+TMEM68


#generate the survival object
data.surv = Surv(finalData$survivalMonth, event=finalData$status)
cox = coxph(data.surv ~ AVPI1+C1QA+C1QB+C3AR1+CCDC116+CD14+COX17+CRLS1+CST7+CTSC+CTSE+CXCL13+DMRT1+DNM2+EIF6+FCER1G+GSTA3+ITGAX+LPHN2+LRPAP1+MCM5+MORN1+NKRF+NRG4+OSMR+PAFAH1B2+PCYOX1+PRAME+RAB11A+RAP2A+RICTOR+SRSF7+TGFBR2+THTPA+TLR2+TMEM117+TMEM68+sex+grade.code+age, data=finalData)
summary(cox)
plot(survfit(cox), main='Multiple regression',xlab='survival (days)',ylab='Survival probability')
#paste summary into excel then use text to columns to split data

#A positive expt(coef) means an increase in one unit of expression means an x increase in hazard (death)
#he stats at the end of the summary are the nul hypothesis that all the coefficients are zero.

#fit a regression linear model to data. Probably invalid due to survival being non-normal
fit = lm(survivalMonth ~ TMEM68, data = finalData)
layout(matrix(c(1,2,3,4),2,2))
plot(fit)