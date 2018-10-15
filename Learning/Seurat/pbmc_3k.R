library(Seurat)
library(dplyr)

#### Load the PBMC dataset ####
pbmc.data <- Read10X(data.dir = "/Volumes/SlowDrive/NGSdata/tenX_genomics/filtered_gene_bc_matrices/hg19/")

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = pbmc.data))
dense.size
sparse.size <- object.size(x = pbmc.data)
sparse.size
dense.size/sparse.size

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
                           project = "10X_PBMC")

#### Perfom the detailed filtering ####

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# GenePlot is typically used to visualize gene-gene relationships, but can be used for anything calculated by the object, 
# there is a rare subset of cells with an outlier level of high mitochondrial percentage and also low UMI we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")

# We filter out cells that have unique gene counts over 2,500 or less than 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate'.  -Inf and Inf should be used if you don't want a lower or upper threshold.
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

#### Normalise the data ####
# Default a global-scaling normalization method “LogNormalize” that normalizes the gene expression measurements 
# for each cell by the total expression, multiplies this by a scale factor (10,000 default), log-transforms the result.
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# Variable genes
par(mfrow = c(1, 1))
# Calculates the average expression and dispersion for each gene, places these genes into bins, 
# and then calculates a z-score for dispersion within each bin. 
# This helps control for the relationship between variability and average expression
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)

#### Removal of technical effects AND PCA ####
# Seurat constructs linear models to predict gene expression based on user-defined variables. 
# The scaled z-scored residuals of these models are stored in the scale.data slot, 
# and used for dimensionality reduction and clustering.
# We could add a cell cycle score but this example doesn't do this.
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = pbmc, pcs.use = 1:3)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)

# Takes a pre-computed PCA (typically calculated on a subset of genes) and projects this onto the entire dataset (all genes). 
# Note that the cell loadings remains unchanged, but now there are gene loading scores for all genes.
pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)

PCHeatmap(object = pbmc, pc.use = 2, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = pbmc, pc.use = 1:12, cells.use = 100, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

#### Determine the significant PCs ####
# To overcome  extensive technical noise in any single gene Seurat clusters cells based on 
# their PCA scores, with each PC essentially representing a ‘metagene’ combines information across a correlated gene set. 
# Determining how many PCs to include downstream is therefore an important step.
pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE)
JackStrawPlot(object = pbmc, PCs = 1:12)
