# Run Toms signature against his breast and colon TCGA data
library(GSVA)
library(gplots)
library(RColorBrewer)
setwd('~/Documents/TomMetabolism/GSEA/input/')

myPalette <- colorRampPalette(c("blue", "white", "yellow"))(n = 1000)


tomsSigs = read.delim('Gene-set_warbug_Literature.gmx', na.strings="")

colon = read.delim('colonCancerPrimary.txt', row.names=1)
colonM = as.matrix(colon)

signatures = list("warburgSet" = tomsSigs$warburgSet, "criticalSet" = tomsSigs$criticalSet, 
                  "DownRegulated" = tomsSigs$Downregulated, "Upregulated"=tomsSigs$Upregulated)

sigs = lapply(signatures, na.omit)

# Run GSVA against TOMs signatures
colonResult = gsva(colonM, signatures,  rnaseq=F, verbose=T, parallel.sz=1)
colonResult = t(colonResult$es.obs)

# Make heat map with Tom's signatures
heatmap.2(t(colonResult), cexRow=1.25, main="Enrichment of manually annotated glycolytic signatures", 
          keysize=1, trace="none", col=myPalette, density.info="none",
          labRow=colnames(colonResult), xlab="Samples", labCol=NA, offsetRow=c(1,1), margins=c(3,9))