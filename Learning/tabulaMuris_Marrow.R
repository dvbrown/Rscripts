require("Matrix")
setwd("/Volumes/SlowDrive/NGSdata/tenX_genomics/tabulaMuris/archive/")

cellbarcodes <- read.table("droplet/Marrow-10X_P7_2/barcodes.tsv")
genenames <- read.table("droplet/Marrow-10X_P7_2//genes.tsv")
molecules <- Matrix::readMM("droplet/Marrow-10X_P7_2//matrix.mtx")

# The cell barcodes are all barcodes and there will be duplicates in other datasets
# Attach the run name
rownames(molecules) <- genenames[,1]
colnames(molecules) <- paste("10X_P7_2", cellbarcodes[,1], sep="_")
meta <- read.delim("droplet_metadata.csv", sep=",", header=TRUE)
head(meta)
# Make formatting consistent
meta[meta$channel == "10X_P7_2",]
mouseID <- "3_56_M"
ann <- read.delim("droplet_annotation.csv", sep=",", header=TRUE)
head(ann)

ann[,1] <- paste(ann[,1], "-1", sep="")
ann_subset <- ann[match(colnames(molecules), ann[,1]),]
celltype <- ann_subset[,3]

# Now lets build the cell-metadata dataframe:
cell_anns <- data.frame(mouse = rep(mouseID, times=ncol(molecules)), type=celltype)
rownames(cell_anns) <- colnames(molecules)

#### The other batch of marrow ####
cellbarcodes2 <- read.table("droplet/Marrow-10X_P7_3//barcodes.tsv")
genenames2 <- read.table("droplet/Marrow-10X_P7_3//genes.tsv")
molecules2 <- Matrix::readMM("droplet/Marrow-10X_P7_3//matrix.mtx")

# The cell barcodes are all barcodes and there will be duplicates in other datasets
# Attach the run name
rownames(molecules2) <- genenames2[,1]
colnames(molecules2) <- paste("10X_P7_3", cellbarcodes2[,1], sep="_")
meta2 <- read.delim("droplet_metadata.csv", sep=",", header=TRUE)
head(meta2)
# Make formatting consistent
meta2[meta2$channel == "10X_P7_3",]
mouseID2 <- "3_57_F"
ann2 <- read.delim("droplet_annotation.csv", sep=",", header=TRUE)
head(ann2)
ann2[,1] <- paste(ann2[,1], "-1", sep="")
ann_subset2 <- ann2[match(colnames(molecules2), ann2[,1]),]
celltype2 <- ann_subset2[,3]
# Now lets build the cell-metadata dataframe:
cell_anns2 <- data.frame(mouse = rep(mouseID2, times=ncol(molecules2)), type=celltype2)
rownames(cell_anns2) <- colnames(molecules2)

#### Sanity checking ####
identical(rownames(molecules), rownames(molecules2))
sum(colnames(molecules) %in% colnames(molecules2))

all_molecules <- cbind(molecules, molecules2)
all_cell_anns <- as.data.frame(rbind(cell_anns, cell_anns2))
all_cell_anns$batch <- rep(c("10X_P7_2", "10X_P7_3"), 
        times = c(nrow(cell_anns), nrow(cell_anns2)))

# Convert to single cell experiment object for disk storage
require("SingleCellExperiment")
require("scater")
all_molecules <- as.matrix(all_molecules)
sceset <- SingleCellExperiment(assays = list(counts = as.matrix(all_molecules)), 
        colData=all_cell_anns)

# Save the object to file
saveRDS(sceset, "/Volumes/SlowDrive/NGSdata/tenX_genomics/tabulaMuris/marrow_droplet.rds")
