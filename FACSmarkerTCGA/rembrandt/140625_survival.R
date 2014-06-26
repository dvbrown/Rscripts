library(survival)
library(coin)
setwd('~/Documents/public-datasets/rembrandt/rembrandt_GBM/processedData/survival/')
source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")

callMarkerSubtype = function (signatureScore, CD133cutoff, CD44cutoff) {
    # Takes a dataframe containing the signature scores and adds a new column that calls FACS marker subtype
    signatureScore$subtype = ""
    signatureScore$subtype = ifelse(signatureScore[,"CD133"] > signatureScore[,"CD44"], "CD133", "CD44")
    # Not having and intermediate case is also better for the Kaplan Myer curve
    #signatureScore$subtype[signatureScore[,"CD133"] < 0 & signatureScore[,"CD44"] < 0] = "doubleNegative"
    signatureScore = sort.dataframe(signatureScore, "subtype")
    signatureScore$subtype = as.factor(signatureScore$subtype)
    return (signatureScore)
}

db <- dbConnect(SQLite(), dbname="~/Documents/public-datasets/rembrandt/rembrandt_GBM/processedData/140624_rembrandtGBM.sqlite")
dbListTables(db)

data = dbReadTable(db, 'facsSubtyeRembrandtProbeMean', row.names='row_names')
#data = dbReadTable(db, 'facsSubtyeRembrandtProbeHighest', row.names='Row_names__1')
clin = dbReadTable(db, 'rembrandtClinical', row.names='patient')

# Try out a double neg subtype
data = callMarkerSubtype(data, 0, 0)

data$verhaak = ""
data$verhaak[data$index == 1] = 'Proneural'
data$verhaak[data$index == 2] = 'Neural'
data$verhaak[data$index == 3] = 'Classical'
data$verhaak[data$index == 4] = 'Mesenchymal'

clinical = clin
clinical$survival = gsub('--', NA,clinical$survival)
clinical$survival = as.numeric(clinical$survival)
clinical$status = 1
clinical$status[is.na(clinical$survival)] = 0

# match up the data and clinical
matched = intersect(row.names(clinical), row.names(data))
dataMatch = data[matched,]
boundData = merge(dataMatch, clinical, by.x='row.names', by.y ='row.names')
boundData$subtype = as.factor(boundData$subtype)
boundData$verhaak = as.factor(boundData$verhaak)
boundData = sort.dataframe(boundData, 9, highFirst=F)

#generate the survival object and plot a Kaplan-Meier
survRembrant = Surv(boundData$survival, event=boundData$status)
surFitRembrandt = survfit(survRembrant~subtype, boundData)

surFitSubtype = survfit(survRembrant~verhaak, boundData)

#### Plot FACS subtype ####
plot(surFitRembrandt, main='FACS marker coexpression signature in Rembrandt Glioblastoma',
     ylab='Survival probability',xlab='survival (months)', 
     col=c("red",'blue'),#'green'),
     cex=1.75, conf.int=F, lwd=1.5)
legend('topright', c('CD133', 'CD44'),# 'Intermediate'), 
       col=c("red",'blue'),#'green'),
       lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

#test for a difference between curves
test = surv_test(survRembrant ~ boundData$subtype)#, subset=!boundData$subtype %in% "intermediate")
test
#legend(locator(1), legend='p = 0.65')

#### Plot verhaak subtype ####
plot(surFitSubtype, main='Verhaak Signature in Rembrandt Glioblastoma',
     ylab='Survival probability',xlab='survival (months)', 
     col=c('green', 'blue', 'orange', 'red'),cex=1.75, conf.int=F, lwd=1.5)

legend('topright', c('Neural', 'Classical', 'Mesenchymal', 'Proneural'), 
       col=c('green', 'blue', 'orange', 'red'), lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

test = surv_test(survRembrant ~ boundData$subtype_y)#, subset=!boundData$subtype %in% "intermediate")
test
#legend(locator(1), legend='p = 0.921')

#### Plot Double Neg ####
surFitDN = survfit(survRembrant~subtype, boundData)

plot(surFitDN, main='FACS marker coexpression signature in Rembrandt Glioblastoma',
     ylab='Survival probability',xlab='survival (months)', 
     col=c("red",'blue','orange'),
     cex=1.75, conf.int=F, lwd=1.5)
legend('topright', c('CD133', 'CD44', 'double negative'), 
       col=c("red",'blue', 'orange'),
       lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

#test for a difference between curves
testDN = surv_test(survRembrant ~ boundData$subtype)#, subset=!boundData$subtype %in% "doubleNegative")
testDN
legend(locator(1), legend='p = 0.34')