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

plotData = read.delim("141107_molecularFeaturesSubtype.txt")
#colors = c("purple", "red", "blue", "green")
colors = c("blue", "red", "green", "purple")

#### Plot molecular subtype per group ####
ggplot(data = plotData[c(1:4, 11:14),], aes(x = subtype, y = percentWithin, fill = entry)) + 
    geom_bar(stat="identity", colour="black") + scale_fill_manual(values=colors) +
    xlab("Coexpression subtype") + ylab("Percent of coexpression subtype") +
    ggtitle("Composition of molecular subtype\nfor coexpression subtype") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))

#### Plot the rest of the data percent between ####
ggplot(data = plotData[c(5:10, 15:20),], aes(x = entry, y = percentBetween, fill = subtype)) + 
    geom_bar(stat="identity", colour="black") + scale_fill_manual(values=c('tomato', 'steelblue1')) +
    xlab("Molecular feature") + ylab("Percent of all patietns") +
    ggtitle("Distribution of molecular features for coexpression subtype") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))

dbDisconnect(db)