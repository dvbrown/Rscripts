#How to load into Seurat ('SuppTable2_PolypALL_merged.txt') and initialize object for analyses presented in Fig. 1 of Ordovas-Montanes et al.,
#NB This analysis was conducted in Seurat v1.4, see here to download package: http://satijalab.org/seurat/get_started_v1_4.html
#Load necessary packages
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(tidyr)
library(dropbead)

#Point to where the files have been downloaded, load in merged cells-by-genes digital gene expression matrix, and inspect top
setwd("~/Projects/Seqwell/Ordovas-Montanes_etal/")
polyp.data=read.table("SuppTable2_PolypALL_merged.txt", header=TRUE, row.names = 1)
head(polyp.data)

plotCumulativeFractionOfReads(polyp.data)

#Initialize seurat object (names.delim set to '_' and Polyp1, Polyp2, etc. refers to biopsy taken, see Supplementary Table 9 for patient characteristics)
polyp <- new("seurat", raw.data = polyp.data)
polyp <- Setup(polyp, min.cells = 5, min.genes = 300, do.logNormalize = T, do.scale= T, do.center = T, names.field = 1, names.delim = "_", project = "Polyp_ALL_TOT")
polyp

#Add metadata to object for polyp status
polyp <- AddMetaData(polyp, NA, col.name = 'polyp')
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp1TOT'] <- 'YES'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp2TOT'] <- 'YES'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp3TOT'] <- 'NO'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp4TOT'] <- 'YES'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp5TOT'] <- 'NO'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp6ATOT'] <- 'NO'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp6BTOT'] <- 'NO'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp7TOT'] <- 'NO'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp8TOT'] <- 'NO'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp9TOT'] <- 'YES'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp11TOT'] <- 'YES'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp12TOT'] <- 'YES'

#Filter cells with UMI counts greater than 12000
polyp <- SubsetData(polyp, subset.name = "nUMI", accept.high = 12000)

#Determine variable genes for input into PCA
polyp <- MeanVarPlot(polyp ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.13, x.high.cutoff = 7, y.cutoff = 0.28, do.contour = T)
length(polyp@var.genes)

#Run PCA
polyp <- PCA(polyp, pc.genes = polyp@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 10)

#Visualize PCA gene loadings and PCA space
VizPCA(polyp,1, num.genes= 60, do.balanced = TRUE)
VizPCA(polyp,2, num.genes= 60, do.balanced = TRUE)
VizPCA(polyp,3, num.genes= 60, do.balanced = TRUE)
VizPCA(polyp,4, num.genes= 60, do.balanced = TRUE)

PCAPlot(polyp, 1, 2)

#Threshold which PCs to use for further dimensionality reduction
PCElbowPlot(polyp)

#Find clusters and run TSNE over the first 12 PCs
polyp <- FindClusters(polyp, pc.use = 1:12, resolution = 1.2, print.output = 0, save.SNN = T)
polyp <- RunTSNE(polyp, dims.use = 1:12, do.fast = T)

#Visualize results by clusters identified, original identity of sample, or polyp status
TSNEPlot(polyp, do.label= T, pt.size = 0.5)
TSNEPlot(polyp, do.label= F, group.by= 'polyp', pt.size = 0.5)
TSNEPlot(polyp, do.label= F, group.by= 'orig.ident', pt.size = 0.2)

#Identify marker genes for clusters using a ROC test and print table
polyp.markers.roc= FindAllMarkers(polyp,return.thresh = 0.65, only.pos = TRUE, test.use = "roc", do.print = TRUE)
write.table(polyp.markers.roc, 'PolypALLTOT_markers_roc.txt', sep = '\t', col.names= NA)

#Display some of the top marker genes for clusters used in Figure 1
FeaturePlot(polyp, 'KRT5', pt.size = 0.7)
FeaturePlot(polyp, 'KRT8', pt.size = 0.7)
FeaturePlot(polyp, 'LTF', pt.size = 0.7)
FeaturePlot(polyp, 'FOXJ1', pt.size = 0.7)
FeaturePlot(polyp, 'TRBC2', pt.size = 0.7)
FeaturePlot(polyp, 'CD79A', pt.size = 0.7)
FeaturePlot(polyp, 'COL1A2', pt.size = 0.7)
FeaturePlot(polyp, 'SPARCL1', pt.size = 0.7)
FeaturePlot(polyp, 'DARC', pt.size = 0.7)
FeaturePlot(polyp, 'HLA-DRA', pt.size = 0.7)
FeaturePlot(polyp, 'TYROBP', pt.size = 0.7)
FeaturePlot(polyp, 'TPSAB1', pt.size = 0.7)
FeaturePlot(polyp, 'CLC', pt.size = 0.7)

#Name cells for Figure 1
polyp <- AddMetaData(polyp, NA, col.name = 'subset')
polyp@data.info$subset[polyp@data.info$res.1.2 == 12] <- 'Basal'
polyp@data.info$subset[polyp@data.info$res.1.2 == 2] <- 'Basal'
polyp@data.info$subset[polyp@data.info$res.1.2 == 8] <- 'Basal'
polyp@data.info$subset[polyp@data.info$res.1.2 == 1] <- 'Apical'
polyp@data.info$subset[polyp@data.info$res.1.2 == 0] <- 'Apical'
polyp@data.info$subset[polyp@data.info$res.1.2 == 4] <- 'Apical'
polyp@data.info$subset[polyp@data.info$res.1.2 == 13] <- 'Glandular'
polyp@data.info$subset[polyp@data.info$res.1.2 == 3] <- 'Glandular'
polyp@data.info$subset[polyp@data.info$res.1.2 == 21] <- 'DOUBLET'
polyp@data.info$subset[polyp@data.info$res.1.2 == 16] <- 'Ciliated'
polyp@data.info$subset[polyp@data.info$res.1.2 == 9] <- 'TCell'
polyp@data.info$subset[polyp@data.info$res.1.2 == 19] <- 'DOUBLET'
polyp@data.info$subset[polyp@data.info$res.1.2 == 20] <- 'DOUBLET'
polyp@data.info$subset[polyp@data.info$res.1.2 == 15] <- 'PlasmaCell'
polyp@data.info$subset[polyp@data.info$res.1.2 == 17] <- 'PlasmaCell'
polyp@data.info$subset[polyp@data.info$res.1.2 == 7] <- 'PlasmaCell'
polyp@data.info$subset[polyp@data.info$res.1.2 == 10] <- 'PlasmaCell'
polyp@data.info$subset[polyp@data.info$res.1.2 == 5] <- 'Fibroblast'
polyp@data.info$subset[polyp@data.info$res.1.2 == 14] <- 'Fibroblast'
polyp@data.info$subset[polyp@data.info$res.1.2 == 6] <- 'Endothelial'
polyp@data.info$subset[polyp@data.info$res.1.2 == 11] <- 'Myeloid'
polyp@data.info$subset[polyp@data.info$res.1.2 == 18] <- 'MastCell'

TSNEPlot(polyp, do.label= F, pt.size = 0.5, group.by = 'subset')

#Removing biological cell doublets detected from analysis of marker gene lists and plot Figure 1B
select0 <- names(polyp@ident[polyp@ident == 0])
select1 <- names(polyp@ident[polyp@ident == 1])
select2 <- names(polyp@ident[polyp@ident == 2])
select3 <- names(polyp@ident[polyp@ident == 3])
select4 <- names(polyp@ident[polyp@ident == 4])
select5 <- names(polyp@ident[polyp@ident == 5])
select6 <- names(polyp@ident[polyp@ident == 6])
select7 <- names(polyp@ident[polyp@ident == 7])
select8 <- names(polyp@ident[polyp@ident == 8])
select9 <- names(polyp@ident[polyp@ident == 9])
select10 <- names(polyp@ident[polyp@ident == 10])
select11 <- names(polyp@ident[polyp@ident == 11])
select12 <- names(polyp@ident[polyp@ident == 12])
select13 <- names(polyp@ident[polyp@ident == 13])
select14 <- names(polyp@ident[polyp@ident == 14])
select15 <- names(polyp@ident[polyp@ident == 15])
select16 <- names(polyp@ident[polyp@ident == 16])
select17 <- names(polyp@ident[polyp@ident == 17])
select18 <- names(polyp@ident[polyp@ident == 18])
select19 <- names(polyp@ident[polyp@ident == 19])
select20 <- names(polyp@ident[polyp@ident == 20])
select21 <- names(polyp@ident[polyp@ident == 21])
nodoublet <- c(select0, select1, select2, select3, select4, select5, select6, select7, select8, select9, select10, select11, select12, select13, select14, select15, select16, select17, select18)
polyp <-SubsetData(polyp, cells.use = nodoublet)
polyp <- SetAllIdent(polyp,'subset') 
polyp@ident <- factor(polyp@ident, levels(polyp@ident)[c(2,1,6,3,5,4,9,10,8,7)])
colors.TSNE <- c("#000066", "#0000FF", "#6699CC", "#99FFFF", '#33CC00', '#00FF66', "#FF9999", '#FF66CC', '#CC00FF', '#CC3333') 
TSNEPlot(polyp, do.label= F, colors.use= colors.TSNE, pt.size= 0.5)

#Toggle identities based on information stored in polyp@data.info, (NB: res.1.2 refers to cluster number identity from FindClusters)
head(polyp@data.info)
polyp <- SetAllIdent(polyp, 'orig.ident')
polyp <- SetAllIdent(polyp, 'res.1.2')

#Subsetting cells and performing differential gene expression of non-polyp vs polyp within a defined subset

#basal epithelial
select <- c(select12, select2, select8)

#basal and secretory and ciliated epithelial
select <- c(select2, select12, select1, select8, select0, select4, select16)

#transitional and secretory epithelial
select <- c(select1, select0, select4)

#secretory epithelial
select <- c(select0, select4)

#glandular epithelial 
select <- c(select3, select13)

#ciliated epithelial 
select <-c(select16)

#Differential expression all cells based on polyp identity
polypepi <-SubsetData(polyp, cells.use = select)
polypepi <- SetAllIdent(polypepi,'polyp') 
select.markers <- FindMarkers(polypepi, 'YES', 'NO', min.pct = 0.1)