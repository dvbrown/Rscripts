library(sqldf)
library(ggplot2)

source("~/Documents/Rscripts/cellBiologyAnalysisFunctions.R")
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

setwd("~/Documents/Cell_biology/141023_summary/")
list.files()

# Intialise and write into database
db <- dbConnect(SQLite(), dbname="assaySummary.sqlite")

# I am going to have to include the background controls if I want to compare across experiments
growth = read.delim("141023_resazurinSummaryRound.txt")
invasion = read.delim("141023_invasionSumRound.txt")
elda = read.delim("141023_eldaSummary.txt")

bw = c("grey21", "grey82", "grey52", "grey97")
colour = c("chartreuse4", "gold", "skyblue2", "orangered1")

################## Growth TMZ assay ############################
# Plot the raw results
growthPlot = ggplot(growth[growth$treatment %in% 'DMSO',], 
                     aes(x=clone, y=mean, fill=cd133status)) + 
    scale_fill_manual(values=colour) +
    #scale_fill_manual(values=bw) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("PDAC") + ylab("Fluorescent intensity") +
    ggtitle("Comparing growth at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))








# Write into database
dbWriteTable(conn = db, name = "resazurin", value = growth, row.names = TRUE)
dbWriteTable(conn = db, name = "invasion", value = invasion, row.names = TRUE)
dbWriteTable(conn = db, name = "elda", value = elda, row.names = TRUE)

dbDisconnect(db)