#a simple Kaplan-Meier
library(survival)
setwd('~/Documents/public-datasets/TCGA/')
clinical = read.delim('./clinicalData/clinical_patient_gbm-sanitised.txt')

extractDataClinical = function(clinicalData) { #A function to get useful clinical data and match rownames to expression
  sur = clinicalData[,c(2,7,8,10,13,19,29)]
  barcode = as.character(clinicalData$bcr_patient_barcode)
  barcode = gsub('$', '.01', barcode)
  sur = cbind(barcode, sur)
  row.names(sur) = sur$barcode
  return (sur)
}

#Arrange data so right censored paitents coded 0 with the days to death changed to days to last followup
censorData = function(survivalData) {
  status = survivalData$vital_status
  status = gsub('LIVING', 0, status)
  status = gsub('DECEASED', 1, status)
  status = as.numeric(status)
  deathFollow = (survivalData[,c(4,5)])
  censored = ifelse(deathFollow$days_to_death %in% NA, deathFollow$days_to_last_followup, deathFollow$days_to_death)
  return (censored)
}

#restructure the survival dataframe
surData = sur[,c(2,6,7)]
surData = cbind(surData, censor, status)
colnames(surData) = c('age.diagnosis', 'gender', 'karnofsky', 'survival', 'death')

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(surData$survival, surData$death)
sur.fit =survfit(data.surv~surData$gender)
plot(sur.fit, main='Kaplan-Meier',ylab='probability',xlab='survival(days)', col=c('red','blue'),xlim=c(0,3000))
legend('topright', c('Female', 'Male'), col=c('blue', 'red'),lwd=1)
summary(data.surv)
#test for a difference between curves
test = survdiff(data.surv~surData$gender)

un.tFrame = t(t.2)
low.creb = subset.data.frame(un.tFrame[,1:91] == range(-1:1))