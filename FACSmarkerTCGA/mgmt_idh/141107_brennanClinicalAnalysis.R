library(sqldf)
library(survival)
library(coin)
library(ggplot2)

#source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")
source("~/Documents/Rscripts/multiplot.R")
db = dbConnect(SQLite(), dbname="~/Documents/public-datasets/cancerBrowser/tcgaData.sqlite")
dbListTables(db)

# Mund and write brennan data into database
setwd("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/mgmtIDHCoexp/")

clinical = dbReadTable(db, "brennanCoexpClinical", row.names=1)
clinical$subtype = as.factor(clinical$subtype)
clinical$Molecular_subtype = as.factor(clinical$Molecular_subtype)
clinical$G_CIMP_status = as.factor(clinical$G_CIMP_status)
clinical$IDH1_status = as.factor(clinical$IDH1_status)
clinical$MGMT_Status = as.factor(clinical$MGMT_Status)

xtabs(~ subtype + Molecular_subtype + G_CIMP_status + IDH1_status + MGMT_Status, data=clinical)


dbDisconnect(db)