getwd()
library(WGCNA)
library(sqldf)
library(ggplot2)

agilent = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt", row.names=1)
rnaseq = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix", row.names=1)
setwd('/Users/d.brown6/Documents/public-datasets/cancerBrowser/deDupAgilent/results/')
db <- dbConnect(SQLite(), dbname='~/Documents/public-datasets/cancerBrowser/deDupAgilent/coexpression.sqlite')

source('~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R')
# The Agilent correlations
cd133Ag = dbReadTable(db, "cd133Allgenes")
cd133Ag$genes = row.names(cd133Ag)
cd44Ag = dbReadTable(db, "cd44Allgenes")
cd44Ag$genes = row.names(cd44Ag)

# Create the RNAseq correlations
cd133Seq = as.data.frame(correlateGeneWithGEM(t(rnaseq), 'PROM1'))
cd133Seq$genes = row.names(cd133Seq)
cd44Seq = as.data.frame(correlateGeneWithGEM(t(rnaseq), 'CD44'))
cd44Seq$genes = row.names(cd44Seq)

cd133 = merge.data.frame(cd133Ag, cd133Seq, by.x='genes', by.y='genes')
cd44 = merge.data.frame(cd44Ag, cd44Seq, by.x='genes', by.y='genes')

# Scatter CD133
ggplot(data=cd133, aes(x=correlation.x, y=correlation.y)) + 
    geom_point(shape=19, alpha=1/4) + geom_smooth(method=lm, colour='red') +
    xlab("Agilent") + ylab("RNAseq") + # Set axis labels
    ggtitle("Correlation of RNAseq and Agilent") +  # Set title
    theme_bw(base_size=18)

# Scatter CD44
ggplot(data=cd44, aes(x=correlation.x, y=correlation.y)) + 
    geom_point(shape=19, alpha=1/4) + geom_smooth(method=lm, colour='red') +
    xlab("Agilent") + ylab("RNAseq") + # Set axis labels
    ggtitle("Correlation of RNAseq and Agilent") +  # Set title
    theme_bw(base_size=18)