# Restart the methylation analysis using cancer browser numeric data
library(sqldf)
library(gplots)
library(RColorBrewer)
setwd("~/Documents/public-datasets/cancerBrowser/methylation/")
list.files()

# Open a connection to database
db <- dbConnect(SQLite(), dbname="~/Documents/public-datasets/cancerBrowser/tcgaData.sqlite")
dbListTables(db)
clinicalRNAseq = dbReadTable(db, "subtypeCallRNAseq")
clinicalAgilent = dbReadTable(db, "subtypeCallAgilent")

# Start with h27k
h27 = read.delim("TCGA_GBM_hMethyl27-2014-08-22/genomicMatrix", row.names =1)
head(h27)
