library(sqldf)
library(ggplot2)

source("~/Documents/Rscripts/cellBiologyAnalysisFunctions.R")
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

setwd("~/Documents/Cell_biology/141023_summary/")
list.files()

# Intialise and write into database
db <- dbConnect(SQLite(), dbname="assaySummary.sqlite")
invasion = read.delim("141023_invasionSumRound.txt")
elda = read.delim("141023_eldaSummary.txt")

bw = c("grey21", "grey82", "grey52", "grey97")
color = c("chartreuse4", "skyblue2", "gold", "orangered1")