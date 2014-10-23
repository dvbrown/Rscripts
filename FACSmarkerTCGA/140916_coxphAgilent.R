library(survival)
library(coin)
library(ggplot2)
library(sqldf)

source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")
source("~/Documents/Rscripts/multiplot.R")

removeEmptyField<- function (dataFrame, columnName) {
  # Remove entries with empty or NA values from dataFrame
  # Column name is a character string
  dataFrame = dataFrame[!is.na(dataFrame[,columnName]),]
  dataFrame = dataFrame[!(dataFrame[,columnName] %in% ""),]
  return (dataFrame)
}

############################ IO ################################
db = dbConnect(SQLite(), dbname="~/Documents/public-datasets/cancerBrowser/tcgaData.sqlite")
dbListTables(db)
clinical = dbReadTable(db, "clinicalAllPatients", row.names=1)
clinical = clinical[!is.na(clinical$X_EVENT),]
row.names(clinical) = clinical$Row_names__1
clinical = removeEmptyField(clinical, "G_CIMP_STATUS")

markerAgilent = dbReadTable(db, "markerScoresAgilent", row.names=1)
row.names(markerAgilent) = gsub("_", ".", row.names(markerAgilent))
molSubtype = c("blue", "red", "green", "purple")

# Extract the clinical data for the Agilent patients
matched = intersect(row.names(clinical), row.names(markerAgilent))
# Subset clinical data for intersect

#################### Do Agilent coxph #########################
clinAgilent = clinical[matched,]
data.surv = Surv(clinAgilent$CDE_survival_time, event=clinAgilent$X_EVENT)
coxPH = coxph(data.surv ~  subtype +  CDE_DxAge  + gender + G_CIMP_STATUS, 
              data=clinAgilent, na.action="na.omit")
summary(coxPH)




dbDisconnect(db)