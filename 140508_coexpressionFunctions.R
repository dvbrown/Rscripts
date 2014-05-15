library(WGCNA)

correlateGeneWithGEM <- function (geneExpressionMatrix = dat, gene='PROM1') {
  # Calculate the correlation between PROM1 expression and all the genes in TCGA GBM
  # The geneExpressionMatrix should be a numeric matrix with genes as columns as rows and patients as rows
  # The gene should be a character string
  corrPval = corAndPvalue(x=geneExpressionMatrix[,gene], y=geneExpressionMatrix)
  
  # Extract the correlation and p-value from the returned list
  correlation = corrPval$cor
  # Measure the  coefficient of determination
  coeffDeter = correlation^2
  pVal = corrPval$p
  fdr = p.adjust(pVal, method='fdr')
  
  result = t(rbind(correlation, coeffDeter, pVal, fdr))
  colnames(result) = c('correlation', ' Rsquared', 'p-value', 'FDR')
  return (result)
}

plotCoexpression <- function (geneCorrelationMatrix, gene) {
  # The geneCorrelationMatrix should be the result of the correlateGeneWithGEM function
  # The gene should be a character string
  par(mfrow=c(2,2))
  hist(geneCorrelationMatrix[,1], main=paste(gene, 'coexpression'), breaks='FD', xlab='Weighted correlation values')
  hist(geneCorrelationMatrix[,2], main=paste(gene, 'coefficient of determination'), breaks='FD', xlab='Weighted correlation values')
  hist(geneCorrelationMatrix[,4], main=paste(gene,'Significance'), breaks='FD', xlab='FDR corrected p-values')
  plot((geneCorrelationMatrix[,2]), -log10(geneCorrelationMatrix[,4]), main=paste(gene, 'Correlation vs significance'), 
       ylab='-log10 FDR p-value', xlab='Absolute correlation', pch=20)
  par(mfrow=c(1,1))
}

makeSquareCoexpressionMatrix <- function (geneCorrelationMatrix, geneExpressionMatrix) {
  # Returns a square adjacency matrix containing the module of genes highly coexpression with the gene of interest
  geneNames = row.names(geneCorrelationMatrix)
  # Build the network adjacency
  # Use the top correlated genes with PROM1 and measure their correlation with the transcriptome
  longMatrix = adjacency(geneExpressionMatrix,  selectCols = geneNames, #for correlation networks only (see below); can be used to select genes whose adjacencies will be calculated. Should be either a numeric vector giving the indices of the genes to be used, or a boolean vector indicating which genes are to be used.
                             type = "unsigned", power = 6, corFnc = "cor", #corOptions = "use = 'p'",
                             distFnc = "dist", distOptions = "method = 'euclidean'")
  
  # Make the adjacency matrix square
  squareMatrix = (longMatrix[colnames(longMatrix),])
  return(squareMatrix)
}

makeDissimilarity <- function (squareMatrix) {
  # Take the square adjacency matrix and return the dissimilarity matrix. This increases the magnitude of the contrasts
  # Calculate the topological overlap matrix
  similarity = TOMsimilarity(squareMatrix, TOMType='unsigned', verbose=3)
  row.names(similarity) = row.names(squareMatrix)
  colnames(similarity) = row.names(squareMatrix)
  # Calculate dissimilarity matrix
  dissTOM = 1-similarity
  return (dissTOM)
}

buildHeatMap <- function (dissimilarityMatrix, gene='PROM1') {
  # There are several methods for branch cutting; our standard method is the Dynamic Tree Cut from the package dynamicTreeCut
  # Module identification using dynamic tree cut. This is the most basic method and returns 3 modules when the cutHeight is 0.999 (default 0.99)
    
  geneTree = flashClust(as.dist(dissimilarityMatrix), method = "average")
  dynamicMods = cutreeDynamic(dendro = geneTree, cutHeight=0.999, method='tree')
  table(dynamicMods)
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  
  # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
  plotTOM = dissimilarityMatrix^6
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA
  title = paste("Network heatmap plot of", gene, "coexpressed genes")
  # Plot the heatmap
  TOMplot(plotTOM, geneTree, dynamicColors, main = title) #, terrainColors=FALSE)
  #,labRow=prom1CgenesNames, ColorsLeft=NA)
  return (dynamicColors)
}

makeMDS <- function (dissimilarityMatrix, moduleColors, gene='CD133') {
  # Make MDS plot using the dissimilarity matrix and the module colors from flash clustering
  par(mfrow=c(1,1))
  cmd1 = cmdscale((dissimilarityMatrix), 3)
  plot(cmd1, col=moduleColors, main = paste('MDS plot of', gene, 'coexpressed genes'), xlab='Most variation', ylab='Second most variation')
}

cytoScapeInput <- function (dissimilarityMatrix, moduleColors) {
# The following R code allow one to specify connection strenghts input to cytoscape
# Select modules based on some measure
modules = c("blue", "brown")
# Select module probes
inModule = is.finite(match(moduleColors, moduleColors))
#modProbes = probes[inModule]
#match1 = match[modProbes, GeneAnnotation$substanceBXH]
modGenes = row.names(dissimilarityMatrix)#[inModule]

# Select the corresponding topological overlap
modTOM = dissimilarityMatrix#[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)

# Export the network into edge and node list files for cytoscape
# cyt = exportNetworkToCytoscape(modTOM, edgeFile=paste("CytoEdge", paste(modules, collapse="-"), ".txt",sep=""),
#                                 nodeFile=paste("CytoNode", paste(modules, collapse="-"), ".txt",sep=""),
#                                weighted=TRUE, threshold=0.02, nodeNames=modGenes, altNodeNames=modGenes,
#                                nodeAttr = moduleColors[inModule])

cyt = exportNetworkToCytoscape(dissimilarityMatrix, edgeFile=paste("CytoEdge", paste(dissimilarityMatrix, collapse="-"), ".txt",sep=""),
                               nodeFile=paste("CytoNode", paste(dissimilarityMatrix, collapse="-"), ".txt",sep=""),
                               weighted=TRUE, threshold=0.02, nodeNames=modGenes, altNodeNames=modGenes,
                               nodeAttr = moduleColors)#[inModule])
return (cyt)
}

cd133_cytoscape = cytoScapeInput(cd133Dissim, cd133Color)
