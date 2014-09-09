library(WGCNA)
source("~/Documents/Rscripts/120704-sortDataFrame.R")

#mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

coVar <- function(x) {
    result = 100*(sd(x) / mean(x))
    return (abs(result))
}

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

buildCorrelationHeatMap <- function (geneExpressionMatrix=dat, adjacencySquareMatrix, gene='PROM1') {
    gene2plot = row.names(adjacencySquareMatrix)
    plotNetworkHeatmap(geneExpressionMatrix, plotGenes=gene2plot, useTOM=FALSE, 
                       networkType="signed", power=1, main=paste("Correlations of ", gene, "mRNA"))
}

buildTOMHeatMap <- function (geneExpressionMatrix=dat, adjacencySquareMatrix, gene='PROM1') {
    gene2plot = row.names(adjacencySquareMatrix)
    plotNetworkHeatmap(geneExpressionMatrix, plotGenes=gene2plot, useTOM=TRUE, 
                       networkType="unsigned", power=4, main=paste("Topology overlap matrix\nof ", gene, "mRNA")) 
}

makeMDS <- function (dissimilarityMatrix, gene='CD133') {
  # Make MDS plot using the dissimilarity matrix and the module colors from flash clustering
  cmd1 = cmdscale((dissimilarityMatrix), 3)
  geneTree = flashClust(as.dist(dissimilarityMatrix), method = "average")
  dynamicMods = cutreeDynamic(dendro = geneTree, cutHeight=0.99, method='tree')
  table(dynamicMods)
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  plot(cmd1, col=dynamicColors, main = paste('MDS plot of', gene, 'coexpressed genes'), xlab='Most variation', ylab='Second most variation')
}

cytoScapeInput <- function (dissimilarityMatrix, #moduleColors, 
                            coexpressedShortList, gene="PROM1") {
# The following R code allow one to specify connection strenghts input to cytoscape.
# CoexpressedShortList is the dataframe which subset the full coexpression data and contains raw correlation values
# Select all module probes

#inModule = is.finite(match(moduleColors, moduleColors))
modGenes = row.names(dissimilarityMatrix)#[inModule]

# Select the corresponding topological overlap
modTOM = dissimilarityMatrix#[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)

# The retreival of gene IDs doesn't work as there are many to one mappings
# modGeneIDs =getBM(filters="hgnc_symbol", 
#             attributes=c("ensembl_gene_id", 'entrezgene',"description", "hgnc_symbol"),
#             values= dimnames(modTOM),
#             mart= mart)

# Export the network into edge and node list files for cytoscape
cyt = exportNetworkToCytoscape(modTOM, edgeFile=paste(gene, "_CytoEdge", ".txt",sep=""),
                                nodeFile=paste(gene, "_CytoNode", ".txt",sep=""),
                               weighted=TRUE, threshold=0.02, nodeNames=modGenes, altNodeNames=NA,#modGeneIDs[,'ensembl_gene_id'],
                               nodeAttr = coexpressedShortList[,'correlation'])
                                   #c(moduleColors[inModule], coexpressedShortList[,'correlation'], coexpressedShortList[,'FDR']))
return (cyt)
}

subsample10times <- function (geneExpressionMatrix=dat, gene="PROM1", iterations=10, statistic="correlation") {
#Subsample the data matrix to validate
  boot_subsample <- function(x, subsample_size) {
# A function to subset the gene expression matrix with replacement. Use the same number of cases
      #' x the data matrix
      #' subsample_size the number of observations (rows) to select at random from x.
      #' a random subsample of the data matrix, x.
         y = x[sample(x = seq_len(nrow(x)), size = subsample_size, replace=TRUE), ]
         return (y)
  }
  
  subSamplePatientsForCorr <- function (expressionMatrix=dat, marker='PROM1') {
# A wrapper for the correlateGeneWithGEM function that returns only the correlation value
    subDat = boot_subsample(geneExpressionMatrix, length(row.names(geneExpressionMatrix)))
    subsamplingCorr = correlateGeneWithGEM(subDat, gene)
    return (subsamplingCorr[,statistic])
  }
# Execute the functions with replacate the number of times specified by iterations  
  result = replicate(iterations, subSamplePatientsForCorr(geneExpressionMatrix, gene), simplify="array")
  return (result)
}

plotResampling = function(resamplingCorrMatrix, resamplingFDRMatrix, originalCoexpressionMatrix, gene="CD133") {
    # Plot resampling metrics
    par(mfrow=c(2,2))
    
    # Plot the distributions of the standard deviations and add the actual data
    hist(apply(resamplingCorrMatrix, 2, sd), breaks='FD', main=paste("Variation in correlation scores \nacross 10 subsamples for", gene), 
        xlab="Standard deviation", col="blue")
    # Add line that signifies real data
    abline(v=sd(originalCoexpressionMatrix[,1]), col='red')

    # Plot the distributions of the FDR and add the actual data
    hist(apply(resamplingFDRMatrix, 2, mean), breaks='FD', main=paste("Variation in FDR scores \nacross 10 subsamples", gene),
        xlab="Mean of FDR", col="forestgreen")
    # Add line that signifies real data
    abline(v=mean(originalCoexpressionMatrix[,4]), col='red')
    
    # Plot the distribution of the data by boxplot and add actual data
    boxColor = c(rep_len("blue", ncol(cd133SubsamplesCorr)), "red")
    boxplot(cbind(resamplingCorrMatrix, originalCoexpressionMatrix[,1]), main=paste("Distribution of correlation scores \nacross 10 subsamples for", gene), 
            col=boxColor, xlab="Subsample", ylab="Correlation")
    par(mfrow=c(1,1))
}

bindGeneExprClinical <- function (clinicalData, subtypedGeneExpression, signatures) {
    # Merges clinical and FACS marker subtyped gene expression information and annotate a color based on Verhaak subtype
    # signatures is a character vector of the signature names
    boundData = merge.data.frame(clinicalData, subtypedGeneExpression, by.x="row.names", by.y="row.names")
    row.names(boundData) = row.names(clinicalData)
    verhaakSubtype = boundData[,c(signatures, "GeneExp_Subtype")]
    verhaakSubtype$colours = "black"
    verhaakSubtype$colours[verhaakSubtype$GeneExp_Subtype == "Proneural"] = "purple"
    verhaakSubtype$colours[verhaakSubtype$GeneExp_Subtype == "Neural"] = "green"
    verhaakSubtype$colours[verhaakSubtype$GeneExp_Subtype == "Classical"] = "blue"
    verhaakSubtype$colours[boundData$GeneExp_Subtype == "Mesenchymal"] = "red"
    return (verhaakSubtype)
}

bindGeneExprCIMPClinical <- function (clinicalData, subtypedGeneExpression, signatures) {
    # Merges clinical and FACS marker subtyped gene expression information and annotate a color based on Verhaak subtype
    # signatures is a character vector of the signature names
    boundData = merge.data.frame(clinicalData, subtypedGeneExpression, by.x="row.names", by.y="row.names")
    row.names(boundData) = row.names(clinicalData)
    verhaakSubtype = boundData[,c(signatures, "GeneExp_Subtype", "G_CIMP_STATUS")]
    verhaakSubtype$colours = "black"
    verhaakSubtype$colours[verhaakSubtype$GeneExp_Subtype == "Proneural"] = "purple"
    verhaakSubtype$colours[verhaakSubtype$GeneExp_Subtype == "Neural"] = "green"
    verhaakSubtype$colours[verhaakSubtype$GeneExp_Subtype == "Classical"] = "blue"
    verhaakSubtype$colours[verhaakSubtype$GeneExp_Subtype == "Mesenchymal"] = "red"
    verhaakSubtype$colours[verhaakSubtype$G_CIMP_STATUS == "G-CIMP"] = "violet"
    # order the heatmap
    verhaakSubtype$order[verhaakSubtype$GeneExp_Subtype == "Proneural"] = 1
    verhaakSubtype$order[verhaakSubtype$GeneExp_Subtype == "Neural"] = 4
    verhaakSubtype$order[verhaakSubtype$GeneExp_Subtype == "Classical"] = 5
    verhaakSubtype$order[verhaakSubtype$GeneExp_Subtype == "Mesenchymal"] = 3
    verhaakSubtype$order[verhaakSubtype$G_CIMP_STATUS == "G-CIMP"] = 2
    return (verhaakSubtype)
}

callMarkerSubtype <- function (signatureScore, CD133cutoff, CD44cutoff) {
    # Takes a dataframe containing the signature scores and adds a new column that calls FACS marker subtype
    signatureScore$subtype = ""
    signatureScore$subtype = ifelse(signatureScore[,"CD133"] > signatureScore[,"CD44"], "CD133", "CD44")
    # Not having and intermediate case is also better for the Kaplan Myer curve
    #signatureScore$subtype[signatureScore[,"CD133"] < -CD133cutoff & signatureScore[,"CD44"] < CD44cutoff] = "intermediate"
    signatureScore = sort.dataframe(signatureScore, "subtype")
    signatureScore$subtype = as.factor(signatureScore$subtype)
    return (signatureScore)
}
