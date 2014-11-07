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

molSubtype = xtabs(~ Molecular_subtype + subtype, data=clinical)
cimp = xtabs (~ G_CIMP_status + subtype, data=clinical)
idh = xtabs(~ IDH1_status + subtype, data=clinical)
mgmt = xtabs(~ MGMT_Status + subtype, data= clinical)

contTable = rbind(molSubtype, cimp, idh, mgmt)
contTable = as.data.frame(contTable[c(1,3:11),])

# Convert to percentage
contTable$cd133Percent = (contTable$CD133 / (contTable$CD133 + contTable$CD44)) * 100
contTable$cd44Percent = (contTable$CD44 / (contTable$CD133 + contTable$CD44)) * 100
write.table(contTable, "./141107_mungThis.txt", sep='\t')

plotData 

ggplot(data = contTable, aes(x = Type.of.Behavior, y = Sample.Size, fill = Stage.of.Change)) + 
    geom_bar()

dbDisconnect(db)