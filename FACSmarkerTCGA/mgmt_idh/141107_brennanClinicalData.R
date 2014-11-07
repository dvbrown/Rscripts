library(survival)
library(coin)
library(ggplot2)
library(sqldf)

source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")
source("~/Documents/Rscripts/multiplot.R")
db = dbConnect(SQLite(), dbname="~/Documents/public-datasets/cancerBrowser/tcgaData.sqlite")
dbListTables(db)

# Mund and write brennan data into database
setwd("~/Documents/public-datasets/TCGA/2013_dataFreeze/")
brennan = read.delim("Cell2013Brennan.txt", row.names=1, na.strings="")
# dbWriteTable(conn = db, name = "brennanClinical", value = brennan, row.names = TRUE)
subtypeAllPatients = dbReadTable(db, "clinicalAllPatients", row.names=T)

subtypeAllPatients$rowName = substr(subtypeAllPatients$Row_names__1, 1, 12)
subtypeName = gsub("\\.", "-", subtypeAllPatients$Row_names__1)
subtypeName1 = substr(subtypeName, 1, 12)

# Remove duplicate patients. These are always 02 ie recurrents
subtypeName3 = subtypeName1[!duplicated(subtypeName1)]
subtypeAllPatients  = subtypeAllPatients [!duplicated(subtypeName1),]
row.names(subtypeAllPatients) = subtypeName3
# Get rid of 1 recurent
subtypeAllPatients = subtypeAllPatients[!subtypeAllPatients$Row_names__1 %in% 'TCGA.06.0171.02',]

 # Merge the data together
bigDF = merge.data.frame(subtypeAllPatients, brennan, by.x=0, by.y=0)
bigDF1 = bigDF[,c("Row_names__1", "CDE_chemo_adjuvant_tmz", "CDE_radiation_any", "subtype", "Secondary_or_recurrent", "Age_at_procedure",
                 "G.CIMP_status", "IDH1_status", "Molecular_subtype", "Therapy_class", "Vital_status", "OS_.days.", "Progression_status", "PFS_.days.")]

# IDH status to mutant
bigDF1$IDH1_status = as.character(bigDF1$IDH1_status)
bigDF1$IDH1_status[bigDF1$IDH1_status %in% c("R132C", "R132G", "R132H")] = 'MT'
bigDF1$IDH1_status = as.factor(bigDF1$IDH1_status)
levels(bigDF1$IDH1_status) = c("WT", "MT")#, "R132C", "R132G", "R132H")