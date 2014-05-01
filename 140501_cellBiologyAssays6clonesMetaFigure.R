setwd('~/Documents/Cell_biology/metaFigure/140501_fig/')

# Just do day 7 for now

growth = read.delim('1400501_day7GrowthMatched.txt')
growth = growth[growth$treatment %in% 'growth',]

tmz = read.delim('140423_day7TMZprocessed.txt')
invasion = read.delim('140501_invasionD7.txt')
elda = read.delim('140501_ELDApercentData.txt')
elda = elda[!elda$Patient %in% c('030a', '034a'),]