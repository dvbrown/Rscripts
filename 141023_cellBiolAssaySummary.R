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
colnames(growth) = c("patient", "assayDate", "subpop", "rep1", "rep2", "rep3", "mean", "sd", "treatment")
growth$sample = paste(growth[,"patient"], growth[,"subpop"], sep="_")
    
invasion = read.delim("141023_invasionSumRound.txt")
elda = read.delim("141023_eldaSummary.txt")

bw = c("grey21", "grey82", "grey52", "grey97")
color = c("chartreuse4", "skyblue2", "gold", "orangered1")

################## Growth TMZ assay ############################
# Plot the raw results
growthPlot = ggplot(growth[growth$treatment %in% 'DMSO',], 
                     aes(x=clone, y=mean, fill=cd133status)) + 
    scale_fill_manual(values=color) +
    #scale_fill_manual(values=bw) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("PDGC") + ylab("Fluorescent intensity") +
    ggtitle("Comparing growth at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
growthPlot


calcGrowthNormalised = function(dataFrame, patientName) {
    # Normalises a patient by the double negative subpopulation set to 1
    # Patient is a charaacter string of patient eg #035
    # Extract all cases of the individual patient
    patient = dataFrame[dataFrame[,"patient"] %in% patientName,]
    # Extract the double negative
    dn = patient[patient[,"subpop"] %in% "CD44-/CD133-",]
    otherSample = patient[!patient[,"subpop"] %in% "CD44-/CD133-",]
#     otherSample$norm1 = otherSample$rep1 / dn$rep1
#     otherSample$norm3 = otherSample$rep3 / dn$rep3
#     otherSample$norm2 = otherSample$rep2 / dn$rep2
    otherSample$normDN = otherSample$mean / dn$mean
    return (otherSample)
}
calcGrowthNormalised(growth, "#020")



# Write into database
dbWriteTable(conn = db, name = "resazurin", value = growth, row.names = TRUE)
dbWriteTable(conn = db, name = "invasion", value = invasion, row.names = TRUE)
dbWriteTable(conn = db, name = "elda", value = elda, row.names = TRUE)

dbDisconnect(db)