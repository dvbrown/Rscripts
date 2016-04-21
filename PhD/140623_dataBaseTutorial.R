# Open a file as a dataframe in R and then write data to a sqllite database
library(sqldf)
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/')

# Using RSQLite's dbConnect() function, a connection to a database is opened. 
# If the named database does not yet exist, one is created. 

# Open a connection to Test.sqlite database
db <- dbConnect(SQLite(), dbname="Test.sqlite")

# Import the data into a database
# dbRemoveTable(db, 'RNAseqBatch1')
dbWriteTable(conn = db, name = "RNAseqBatch1", value = df, row.names = TRUE)

dbListTables(db)                 # The tables in the database
dbListFields(db, "RNAseqBatch1")       # The columns in a table
head(dbReadTable(db, "RNAseqBatch1"))        # The data in a table

#The connection to the database is closed, and as a precaution, the data frame is removed from R's environment.
dbDisconnect(db)            # Close connection
rm(df, df1)

################### Read in a file directly and exprot to a database ####
db1 <- dbConnect(SQLite(), dbname="Test1.sqlite")

read.csv.sql('GLMedgeR/131021_normalisedCPM.txt', sql="CREATE TABLE normCPM AS SELECT * FROM file",
             dbname="Test1.sqlite", sep='\t')
dbDisconnect(db1) 