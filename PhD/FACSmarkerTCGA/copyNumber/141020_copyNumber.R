library(sqldf)

setwd("~/Documents/public-datasets/TCGA/copyNumber/")
list.files()

############################ IO ################################
anno = read.delim("TCGA_GBM_PCT_CNA_Annotation.txt", row.names=1)
triplet = read.delim("TCGA_GBM_PCT_CNA_Physical_triplet.txt")
summary = read.delim("TCGA_GBM_PCT_CNA_Summary.txt")

db <- dbConnect(SQLite(), dbname="~/Documents/public-datasets/cancerBrowser/tcgaData.sqlite")
dbListTables(db)
clinical = dbReadTable(db, "clinicalAllPatients", row.names=1)
clinical$patient = substr(row.names(clinical), 1, 12)

############################ Separate out the CD133 and CD44 subtypes and take the column sum of each ################################
cd133Patients = clinical[clinical$subtype %in% "CD133",]$patient
cd44Patients = clinical[clinical$subtype %in% "CD44",]$patient
cd133Copy = summary[,cd133Patients]


head(summary)
p = colnames(summary)
mungData = merge.data.frame(summary, clinical, by.x='row.names', by.y="patient")

# Draw scatterplot
ggplot(data=signatureScores, aes(x=CD133, y=CD44, color=Patient)) + 
    geom_point(shape=19, alpha=1) + geom_smooth(method=lm, colour='black') +
    scale_colour_manual(values=cbPalette) +
    xlab("CD133 signatures") + ylab("CD44 signatures") + # Set axis labels
    ggtitle("Anoop et al 2014 single cell RNAseq\nall patients by coexpression signature score") +  # Set title
    theme_bw(base_size=18)