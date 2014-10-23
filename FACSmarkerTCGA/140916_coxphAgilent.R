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
clin = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)
clin = clin[, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender", 'CDE_chemo_adjuvant_tmz', 'CDE_chemo_tmz',
                           'CDE_radiation_any', 'CDE_tmz_chemoradiation_standard', 'GeneExp_Subtype')]

clinical = clinical[!clinical$GeneExp_Subtype %in% "",]

markerAgilent = dbReadTable(db, "markerScoresAgilent", row.names=1)
molSubtype = c("blue", "red", "green", "purple")

# Extract the clinical data for the Agilent patients
matched = intersect(row.names(clinical), row.names(markerAgilent))
# Subset clinical data for intersect


#################### bind the clinical and subtyping info together #########################