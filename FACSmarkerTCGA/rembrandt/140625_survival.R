library(survival)
library(coin)
setwd('~/Documents/public-datasets/rembrandt/rembrandt_GBM/processedData/survival/')

callMarkerSubtype = function (signatureScore, CD133cutoff, CD44cutoff) {
    # Takes a dataframe containing the signature scores and adds a new column that calls FACS marker subtype
    signatureScore$subtype = ""
    signatureScore$subtype = ifelse(signatureScore[,"CD133"] > signatureScore[,"CD44"], "CD133", "CD44")
    # Not having and intermediate case is also better for the Kaplan Myer curve
    signatureScore$subtype[signatureScore[,"CD133"] < 0 & signatureScore[,"CD44"] < 0] = "doubleNegative"
    signatureScore = sort.dataframe(signatureScore, "subtype")
    signatureScore$subtype = as.factor(signatureScore$subtype)
    return (signatureScore)
}

db <- dbConnect(SQLite(), dbname="~/Documents/public-datasets/rembrandt/rembrandt_GBM/processedData/140624_rembrandtGBM.sqlite")
dbListTables(db)

data = dbReadTable(db, 'facsSubtyeRembrandtProbeMean', row.names='Row_names__1')
#data = dbReadTable(db, 'facsSubtyeRembrandtProbeHighest', row.names='Row_names__1')
clin = dbReadTable(db, 'rembrandtClinical', row.names='patient')

# Try out a double neg subtype
data = callMarkerSubtype(data, 0, 0)

# Mung the clinical data for observations
clinical = clin
clinical$survival = gsub('--', NA,clinical$survival)
clinical$survival = as.numeric(clinical$survival)
clinical$status = 1
clinical$status[is.na(clinical$survival)] = 0

# match up the data and clinical
matched = intersect(row.names(clinical), row.names(data))
dataMatch = data[matched,]
boundData = merge(dataMatch, clinical, by.x='row.names', by.y ='row.names')
boundData$subtype_x = as.factor(boundData$subtype_x)
boundData$subtype_y = as.factor(boundData$subtype_y)

#generate the survival object and plot a Kaplan-Meier
survRembrant = Surv(boundData$survival, event=boundData$status)
surFitRembrandt = survfit(survRembrant~subtype_x, boundData)

surFitSubtype = survfit(survRembrant~subtype_y, boundData)

#### Plot FACS subtype ####
plot(surFitRembrandt, main='FACS marker coexpression signature in Rembrandt Glioblastoma',
     ylab='Survival probability',xlab='survival (months)', 
     col=c("red",'blue'),#'green'),
     cex=1.75, conf.int=F, lwd=1.5)
legend('topright', c('CD133', 'CD44'),# 'Intermediate'), 
       col=c("red",'blue'),#'green'),
       lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

#test for a difference between curves
test = surv_test(survRembrant ~ boundData$subtype_x)#, subset=!boundData$subtype %in% "intermediate")
test
#legend(locator(1), legend='p = 0.65')

#### Plot verhaak subtype ####
plot(surFitSubtype, main='Verhaak Signature in Rembrandt Glioblastoma',
     ylab='Survival probability',xlab='survival (months)', 
     col=levels(boundData$subtype_y),cex=1.75, conf.int=F, lwd=1.5)

legend('topright', c('Classical', 'Neural', 'Mesenchymal', 'Pronerual'), col=levels(boundData$subtype_y),
lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

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