library(survival)
setwd('~/Documents/stemCellSig/130117_signature/')
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
  colnames(data) = c(sig, 'barcode', 'patient','age', 'sex', 'KF.score','survivalDay','status')
  data$survivalDay = as.numeric(data$survivalDay) 
  return(data)
}

#prepare the data for multiple regression
data = read.delim('~/Documents/public-datasets/TCGA/expressionData/Expression-Genes/121109_AgilentPooledzTransform.txt', sep='\t',row.names=1)
#stem cell signature change this!
signature = read.delim('~/Documents/stemCellSig/130117_signature/130117_stemCellSig.txt')
signature = as.character(signature[,1]) #the score column ranks the number of sites in the given promoter

clinical = read.delim('~/Documents/public-datasets/TCGA/clinicalData/120731-survivalDataStructered.txt')
clinical = extractDataClinical(clinical)
sig = intersect(row.names(data),signature)
subsetGrav = data[sig,]
data = t(subsetGrav)
row.names(data) = gsub('$', '.CEL', row.names(data))
finalData = bindGEMsurvival(data, clinical)

#got to comment out the geneList for reuse
#LembergerOnlyTop50 = AVPI1+C1QA+C1QB+C3AR1+CCDC116+CD14+COX17+CRLS1+CST7+CTSC+CTSE+CXCL13+DMRT1+DNM2+EIF6+FCER1G+GSTA3+ITGAX+LPHN2+LRPAP1+MCM5+MORN1+NKRF+NRG4+OSMR+PAFAH1B2+PCYOX1+PRAME+RAB11A+RAP2A+RICTOR+SRSF7+TGFBR2+THTPA+TLR2+TMEM117+TMEM68
#FOXG1+GFAP+ID1+NANOG+NES+NOTCH1+OLIG2+PAX6+POU3F2+POU5F1+PROM1+SOX2

#generate the survival object
data.surv = Surv(finalData$survivalDay, event=finalData$status)
cox = coxph(data.surv ~ FOXG1+GFAP+ID1+NANOG+NES+NOTCH1+OLIG2+PAX6+POU3F2+POU5F1+PROM1+SOX2+age+sex+KF.score, data=finalData)
summary(cox)
plot(survfit(cox), main='Multiple regression',xlab='survival (days)',ylab='Survival probability')
#paste summary into excel then use text to columns to split data

#A positive expt(coef) means an increase in one unit of expression means an x increase in hazard (death)
#he stats at the end of the summary are the nul hypothesis that the coefficients are zero.