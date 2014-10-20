# Examine the relationship between G-CIMP and CD133 subtype as well as absolute expression
library(sqldf)
library(ggplot2)

setwd("~/Documents/public-datasets/cancerBrowser/cd133_gCIMP/")

############################ IO ################################
db = dbConnect(SQLite(), dbname="~/Documents/public-datasets/cancerBrowser/tcgaData.sqlite")
dbListTables(db)
clinical = dbReadTable(db, "clinicalAllPatients", row.names=1)
agilentGem = dbReadTable(db, "AgilentGem", row.names=1)

rnaSeqGem = dbReadTable(db, "RNAseqGem", row.names=1)
markerScore = dbReadTable(db, "markerScoresRNAseq", row.names=1)

############################ Annotate the CIMPs ###############################

markerCIMP = merge.data.frame(clinical, markerScore, by.x=0, by.y=0)
markerCIMP = markerCIMP[!markerCIMP$G_CIMP_STATUS %in% "",]

ggplot(markerCIMP, aes(x=G_CIMP_STATUS, y=CD133, fill=G_CIMP_STATUS)) + geom_boxplot() +
    scale_colour_manual(values=c("aqua", "orange")) +
    xlab("CD133 signatures") + ylab("CD44 signatures") + # Set axis labels
    ggtitle("No difference in CD133 coexpression signature by CIMP status") +  # Set title
    guides(fill=FALSE) + geom_jitter() +  theme_bw(base_size=18)

############################ Extract the raw expression values ###############################
markermRNA = rnaSeqGem[,c('PROM1', 'FUT4', 'CD44')]
gemCIMP = merge.data.frame(clinical, markermRNA, by.x=0, by.y=0)
gemCIMP = gemCIMP[!gemCIMP$G_CIMP_STATUS %in% "",]

ggplot(gemCIMP, aes(x=G_CIMP_STATUS, y=PROM1, fill=G_CIMP_STATUS)) + geom_boxplot() +
    scale_colour_manual(values=c("aqua", "orange")) +
    xlab("CD133 signatures") + ylab("CD44 signatures") + # Set axis labels
    ggtitle("Minor difference in CD133 coexpression signature by CIMP status") +  # Set title
    guides(fill=FALSE) + geom_jitter() +  theme_bw(base_size=18)