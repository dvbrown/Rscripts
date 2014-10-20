# Examine the relationship between G-CIMP and CD133 subtype as well as absolute expression
library(sqldf)

setwd("~/Documents/public-datasets/cancerBrowser/cd133_gCIMP/")

############################ IO ################################
db = dbConnect(SQLite(), dbname="~/Documents/public-datasets/cancerBrowser/tcgaData.sqlite")
dbListTables(db)
clinical = dbReadTable(db, "clinicalAllPatients", row.names=1)
agilentGem = dbReadTable(db, "AgilentGem", row.names=1)

rnaSeqGem = dbReadTable(db, "RNAseqGem", row.names=1)
markerScore = dbReadTable(db, "markerScoresRNAseq", row.names=1)

############################ Subset the CIMPs ################################
gcimp = row.names(clinical[clinical$G_CIMP_STATUS %in% "G-CIMP", ])
notGcimp = row.names(clinical[!clinical$G_CIMP_STATUS %in% "G-CIMP", ])

markerCIMP = markerScore[gcimp,]
markerNotcIMP = markerScore[notGcimp,]