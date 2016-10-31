setwd("Dropbox/Manuscripts/PlosONE_majorAmendments/Figures/Figure2/mouse/")
list.files()

dat = read.table("161029_immunoMouseCounts.txt", row.names = 1, sep="\t", header = T)
dat