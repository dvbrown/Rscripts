library(DropletUtils)
library(scater)
library(scran)
library(EnsDb.Hsapiens.v86)

fname <- "/Volumes/SlowDrive/NGSdata/tenX_genomics/filtered_gene_bc_matrices/hg19/"
sce <- read10xCounts(fname, col.names=TRUE)
sce

# Change the gene symbols to gene IDs
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
head(rownames(sce))

# Get the chromosome from where the gene comes from
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, 
                   column="SEQNAME", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="MT")

#### Separate the cell drops from the empty drops ####
bcrank <- barcodeRanks(counts(sce))
# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=bcrank$inflection, col="darkgreen", lty=2)
abline(h=bcrank$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# Set the empty drop threshold to 100 then perform a statistical test (FDR 0.01) to call cell drops from the empties
set.seed(100) # Monte Carlo test is used therefore set seed
e.out <- emptyDrops(counts(sce), lower = 100)
sum(e.out$FDR <= 0.01, na.rm=TRUE)
sum(e.out$FDR >= 0.01, na.rm=TRUE)
# The number of empty drops is printed

table(Sig=e.out$FDR <= 0.01, Limited=e.out$Limited)
# Limited: = Logical, indicating whether a lower p-value could be obtained by increasing npts.

# Subset the single cell experiment so it only contains cells
sce <- sce[,which(e.out$FDR <= 0.01)]

#### Quality control on the cells ####
# This is a very relaxed filter, Aaron recommends going back and forth from clustering to make QC tighter
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(location=="MT")))
par(mfrow=c(1,3))
hist(sce$log10_total_counts, breaks=20, col="grey80",
     xlab="Log-total UMI count")
hist(sce$log10_total_features_by_counts, breaks=20, col="grey80",
     xlab="Log-total number of expressed features")
hist(sce$pct_counts_Mito, breaks=20, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")

# Mitochondrial genes
high.mito <- isOutlier(sce$pct_counts_Mito, nmads=3, type="higher") # nmads means median abs deviation
sce <- sce[,!high.mito]
summary(high.mito)

#### Looking at gene expression ####
ave <- calcAverage(sce)
rowData(sce)$AveCount <- ave
par(mfrow=c(1,1))
hist(log10(ave), col="grey80")

#plotHighestExprs(sce)

#### Normalizing for cell-specific biases ####
# Perform preclustering before normalisation to avoid grouping cells that are too different together in the normalisation
clusters <- quickCluster(sce, method="igraph", min.mean=0.1,
        irlba.args=list(maxit=1000)) # for convergence.
table(clusters)
sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
summary(sizeFactors(sce))

# Look at the correlation between library size and the scaling factor used in normalisation
plot(sce$total_counts, sizeFactors(sce), log="xy")
sce <- normalize(sce)

#### Modelling the mean-variance trend ####
# Fit a trend on technical noise using a poisson distribution
new.trend <- makeTechTrend(x=sce)

fit <- trendVar(sce, use.spikes=FALSE, loess.args=list(span=0.05))
plot(fit$mean, fit$var, pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)
curve(new.trend(x), col="red", add=TRUE)

# We decompose the variance for each gene using the Poisson-based trend
# examine the genes with the highest biological components.
fit0 <- fit
fit$trend <- new.trend
dec <- decomposeVar(fit=fit)
top.dec <- dec[order(dec$bio, decreasing=TRUE),] 
head(top.dec)

# Plot genes with highest biological variation
plotExpression(sce, features=rownames(top.dec)[1:10])

# Dimension reduction
# Decide on how many dimensions to keep of the PCA
sce <- denoisePCA(sce, technical=new.trend, approx=TRUE)
ncol(reducedDim(sce, "PCA"))

plot(attr(reducedDim(sce), "percentVar"), xlab="PC",
     ylab="Proportion of variance explained")
abline(v=ncol(reducedDim(sce, "PCA")), lty=2, col="red")

plotPCA(sce, ncomponents=3, colour_by="log10_total_features_by_counts")

# DO a t-SNE on the PCA reduced spaceq()
sce <- runTSNE(sce, use_dimred="PCA", perplexity=30, rand_seed=100)
plotTSNE(sce, colour_by="log10_total_features_by_counts")

# Clustering with graph-based methods
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)
sce$Cluster <- factor(clusters$membership)
table(sce$Cluster)

cluster.mod <- clusterModularity(snn.gr, sce$Cluster, get.values=TRUE)
log.ratio <- log2(cluster.mod$observed/cluster.mod$expected + 1)
# Modukarity tests if there are random edges in the graph

library(pheatmap)
pheatmap(log.ratio, cluster_rows=FALSE, cluster_cols=FALSE, 
         color=colorRampPalette(c("white", "blue"))(100))
plotTSNE(sce, colour_by="Cluster")

#### Find markers ####
markers <- findMarkers(sce, clusters=sce$Cluster, direction="up")
marker.set <- markers[["2"]]
head(marker.set[,1:8], 10) # only first 8 columns, for brevity

chosen <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=chosen, exprs_values="logcounts", 
            zlim=5, center=TRUE, symmetric=TRUE, cluster_cols=FALSE,
            colour_columns_by="Cluster", columns=order(sce$Cluster))