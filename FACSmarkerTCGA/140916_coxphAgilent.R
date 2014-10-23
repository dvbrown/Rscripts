library(survival)
library(coin)
library(ggplot2)
library(sqldf)

source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")
source("~/Documents/Rscripts/multiplot.R")

############################ IO ################################
db = dbConnect(SQLite(), dbname="~/Documents/public-datasets/cancerBrowser/tcgaData.sqlite")
dbListTables(db)
clinical = dbReadTable(db, "clinicalAllPatients", row.names=1)
clinical = clinical[!is.na(clinical$X_EVENT),]

markerAgilent = dbReadTable(db, "markerScoresAgilent", row.names=1)
molSubtype = c("blue", "red", "green", "purple")

# Extract the clinical data for the Agilent patients
matched = intersect(row.names(clinical), row.names(markerAgilent))
# Subset clinical data for intersect

#################### bind the clinical and subtyping info together #########################


dbDisconnect(db)