#See what overlap just the expression CD133, CD44 and CD133 alone have with Verhaak subtype
library(sqldf)

getwd()
source('/Users/d.brown6/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R')
setwd('~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/signatureComparison/useFACSgeneOnly/')
list.files()

db <- dbConnect(SQLite(), dbname="~/Documents/public-datasets/cancerBrowser/tcgaData.sqlite")

dbListTables(db)                 # The tables in the database
dbListFields(db, "AgilentGem")       # The columns in a table