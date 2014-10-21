# Examine the relationship between G-CIMP and CD133 subtype as well as absolute expression
library(sqldf)
library(ggplot2)
source("~/Documents/Rscripts/multiplot.R")
setwd("~/Documents/public-datasets/cancerBrowser/cd133_gCIMP/")

############################ IO ################################
db = dbConnect(SQLite(), dbname="~/Documents/public-datasets/cancerBrowser/tcgaData.sqlite")
dbListTables(db)
clinical = dbReadTable(db, "clinicalAllPatients", row.names=1)
clinical = clinical[!clinical$GeneExp_Subtype %in% "",]
rnaSeqGem = dbReadTable(db, "RNAseqGem", row.names=1)
markerScore = dbReadTable(db, "markerScoresRNAseq", row.names=1)
molSubtype = c("blue", "red", "green", "purple")

############################ Annotate the CIMPs ###############################
markerCIMP = merge.data.frame(clinical, markerScore, by.x=0, by.y=0)
markerCIMP = markerCIMP[!markerCIMP$G_CIMP_STATUS %in% "",]
sig = ggplot(markerCIMP, aes(x=G_CIMP_STATUS, y=CD133, fill=G_CIMP_STATUS)) + geom_boxplot() +
    scale_colour_manual(values=molSubtype) +
    xlab("RNAseq") + ylab("enrichment score") + # Set axis labels
    ggtitle("No difference in CD133 coexpression\n signature by CIMP status") + # Set title
    guides(fill=FALSE) + geom_jitter(aes(colour=GeneExp_Subtype, shape=GeneExp_Subtype)) + theme_bw(base_size=18)

############################ Extract the raw expression values ###############################
markermRNA = rnaSeqGem[,c('PROM1', 'FUT4', 'CD44')]
gemCIMP = merge.data.frame(clinical, markermRNA, by.x=0, by.y=0)
gemCIMP = gemCIMP[!gemCIMP$G_CIMP_STATUS %in% "",]
mRNA = ggplot(gemCIMP, aes(x=G_CIMP_STATUS, y=PROM1, fill=G_CIMP_STATUS)) + geom_boxplot() +
    scale_colour_manual(values=molSubtype) +
    xlab("RNAseq") + ylab("mRNA expression") + # Set axis labels
    ggtitle("Minor difference in CD133 mRNA expression\n by CIMP status") + # Set title
    guides(fill=FALSE) + geom_jitter(aes(colour=GeneExp_Subtype, shape=GeneExp_Subtype)) + theme_bw(base_size=18)
#multiplot(sig, mRNA, cols=2)

############################ Repeat with Agilent data ###############################
agilentGem = dbReadTable(db, "AgilentGem", row.names=1)
agilentScore = dbReadTable(db, "markerScoresAgilent", row.names=1)
row.names(agilentScore) = gsub("_", ".", row.names(agilentScore))
colnames(agilentGem) = gsub("_", ".", colnames(agilentGem))
markerCIMP = merge.data.frame(clinical, agilentScore, by.x=0, by.y=0)
markerCIMP = markerCIMP[!markerCIMP$G_CIMP_STATUS %in% "",]
markermRNA = t(agilentGem[c('PROM1', 'FUT4', 'CD44'),])
gemCIMP = merge.data.frame(clinical, markermRNA, by.x=0, by.y=0)
gemCIMP = gemCIMP[!gemCIMP$G_CIMP_STATUS %in% "",]

sigAgilent = ggplot(markerCIMP, aes(x=G_CIMP_STATUS, y=CD133, fill=G_CIMP_STATUS)) + geom_boxplot() +
    scale_colour_manual(values=molSubtype) +
    xlab("Agilent") + ylab("enrichment score") + # Set axis labels
    ggtitle("No difference in CD133 coexpression\n signature by CIMP status") + # Set title
    guides(fill=FALSE) + geom_jitter(aes(colour=GeneExp_Subtype, shape=GeneExp_Subtype)) + theme_bw(base_size=18)

mRNAAgilent = ggplot(gemCIMP, aes(x=G_CIMP_STATUS, y=PROM1, fill=G_CIMP_STATUS)) + geom_boxplot() +
    scale_colour_manual(values=molSubtype) +
    xlab("Agilent") + ylab("mRNA expression") + # Set axis labels
    ggtitle("Difference in CD133 mRNA expression\n by CIMP status") + # Set title
    guides(fill=FALSE) + geom_jitter(aes(colour=GeneExp_Subtype, shape=GeneExp_Subtype)) + theme_bw(base_size=18)

multiplot(sigAgilent, mRNAAgilent, cols=2)

dbDisconnect(db)