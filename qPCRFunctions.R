library(reshape)
library(ggplot2)
library(plyr)

# Intialise the package at the end by building a list containing all the functions in this script

# Transpose the 384 well map (tab text) using a python script and output into another file
transposeLinear = function(well384Map, linearMapFile='output.txt') {
  pythonCall = paste('~/Documents/Eclipseworkspace/Bioinformatics/Filtering/transposeLinear.py -i', well384Map, '>', linearMapFile, sep=' ')
  system(pythonCall)
  sampleLabels = read.delim(linearMapFile, header=F)
  colnames(sampleLabels) = c('location', 'sample')
  return (sampleLabels)
}

#split sample names by whitespace. Takes a dataFrame and splits anything with whitespace into 2 columns
splitSampleName = function(plateMap) {
  # The column with sample is the vector containing the sample names you wish to split
  splitted = colsplit.factor(plateMap[['sample']], split = " ", names = c('origin', 'gene'))
  result = cbind(plateMap, splitted)
  return (result)
}

buildDataFrameFromddCT = function(plateMap, CtData) {
  #Bind the dataframes containing the sample labels with the raw data itself and remove useless columns
  rawData = CtData[,c(3,4,5,8)]
  result = merge.data.frame(plateMap, rawData, by.x='location', by.y='Pos')
  return (result)
}

extractReplicates <- function (indexes, ctData) {
  # Retreive the indexes of the 384 wellplate
  indexes = c(1:384)
  
  # Keep only cases with data in them as the merge function doesn't work with NAs
  CtData = ctData[complete.cases(ctData[,3]),]
  # Subset each Cp into its replicates. Takes a vector with the indexes to to subset and then takes the
  # even entries and odd entries separately from the dataframe containing cp values
  even = as.character(indexes[indexes%%2 == 0])
  odd = as.character(indexes[indexes%%2 == 1])
  even = paste('Sample', even, sep=' ')
  odd = paste('Sample', odd, sep=' ')
  
  #rep1 = CtData[odd, c(1:6)]
  rep1 = CtData[CtData$Name %in% odd, c(1:6)]
  
  #rep1 = rep1[complete.cases(rep1$sample),]
  
  #rep2 = CtData[even, c(1:6)]
  rep2 = CtData[CtData$Name %in% even, c(1:6)]
  
  #rep2 = rep2[complete.cases(rep2$sample),]
  
  boundData = merge(rep1, rep2, by.x='sample', by.y='sample')
  ################ Remove columns that do not add information
  usefulData = boundData[,c(1,2,3,4,6,7,11)]
  # Compute the mean and the standard deviation of the replicates
  usefulData$meanCP = rowMeans(cbind(usefulData$Cp.x, usefulData$Cp.y), na.rm=T)
  usefulData$stdDevCP = apply(cbind(usefulData$Cp.x, usefulData$Cp.y), 1, sd, na.rm=T)
  # Package the output in a list
  result = list(rep1, rep2, usefulData)
  return (result)
}

ddCTcalculate = function(geneOfInterest, sampleOfInterest='020_N', houseKeepingGene='GAPDH', referenceSample='020_N', data=rawData) {
  # The gene of interest and the sample of interest are both vectors from the dataFrame.
  # House keeping gene and 
  
  # This function has been successfully debugged and has been shown to work.
  sampleHouse = paste(sampleOfInterest, houseKeepingGene)
  SampleGene = paste(sampleOfInterest, geneOfInterest)
  # Extract the Cp of the hosue keeping gene
  houseCp = data[sampleHouse, 'meanCP']
  # place holder 020_N ATP5G3 as the current row
  geneCp = data[SampleGene, 'meanCP']
  # dCt calculation for the sample of interest
  dCt = houseCp - geneCp
  
  # Extract the meanCP for the reference sample. First get the index of the housekeeping gene, then the gene of interest
  refDctRowHouse = paste(referenceSample, houseKeepingGene)
  refDctRowGene = paste(referenceSample, geneOfInterest)
  # Calculate dCt for the reference sample
  referenceSample_dCt = data[refDctRowHouse, 'meanCP'] - data[refDctRowGene, 'meanCP']
  
  # Calculate ddCt
  ddCtNotSquared = dCt - referenceSample_dCt
  ddCt = 2^ddCtNotSquared
  
  return (ddCt)
}

plot_ddCt = function(Expressionformula, dataFrame, graphTitle='A grouped barchart', yaxisLabel='y axis') {
  # This will make barcharts without error bars
  # Expression formula is of the type ddCt ~ cell type or whatever you want the bars to be grouped by
  p = barchart(Expressionformula, data = dataFrame, groups = gene.x, 
                scales = list(x = list(rot=90,cex=0.8)), main = graphTitle, ylab=yaxisLabel,
                auto.key=list(space="top", columns=3,
                title="genes", cex.title=1))      
  #returns a plot object that when you look at it plots stuff
  return (p)
}

niceGroupedBarPlot <- function (dataFrame, ddCt='ddCt', sampleOrigin="origin.x", gene="gene.x", graphTitle="A pretty plot") {
  ggplot(data=dataFrame, aes(x=sampleOrigin, y=ddCt, fill=gene)) + 
      geom_bar(stat="identity", position=position_dodge(), colour="black") + 
      scale_fill_hue(name="Gene") +      # Set legend title
      xlab("Sample") + ylab("ddCt") + # Set axis labels
      ggtitle(graphTitle) +  # Set title
      theme_bw(base_size=18)
}

niceErrorBarPlot <- function (summarisedData, xAxis=gene.x, yAxis=mean, groupVariable=cd133, 
                              title='A title', xLabel='Gene', yLabel='Expression') {
  p = ggplot(summarisedData, aes(x=xAxis, y=yAxis, fill=groupVariable)) + 
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=yAxis-se, ymax=yAxis+se),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9)) +
        xlab(xLabel) +
        ylab(yLabel) +
        scale_fill_hue(name="CD133")+#, Legend label, use darker colors
        ggtitle(title) +
        scale_y_continuous(breaks=0:20*4) +
        # Setting vjust to a negative number moves the asterix up a little bit to make the graph prettier
        geom_text(aes(label=star), colour="black", vjust=-2, size=10) +
        theme_bw(base_size=16)
  return (p)
}

build_ddCTmatrix = function(ddCtFile, originColumn=2, geneColumn=3, ddCtColumn=9, output='matrix.txt') {
    # A function to coerce ddCT values and genes into a double matrix for statistical analysis
    # origin is the name of the sample eg #020 CD133 negative
    pythonCall = paste('./buildNumericMatrix_ddCt.py', '-i', ddCtFile,
                       '-o', originColumn, '-g', geneColumn, '-d', ddCtColumn,
                       '>', output, sep=' ')
    system(pythonCall)

    f = read.delim('matrix.txt', header=T)
    # sort the dataframe
    g = f[with(f, order(f[,1], decreasing=F)),]
    # set rownames then remove
    row.names(g) = g[,1]
    return (g)
}

summariseStatistics_ddCt <- function (dataFrame, groupVariableA='cd133', gropVariableB='gene.x') {
    # Generate N, mean, sd and se statistics for a dataframe
    cData = ddply(cd133negPos, c(groupVariableA, gropVariableB), summarise,
                N    = sum(!is.na(ddCt)),
                mean = mean(ddCt, na.rm=TRUE),
                sd   = sd(ddCt, na.rm=TRUE),
                se   = sd / sqrt(N) )
    
    # Add a column with stars describing if a test is significant
#     percentData$star <- " "
#     percentData$star[percentData$adjust < .05]  = "*"
#     percentData$star[percentData$adjust < .01]  <- "**"
#     percentData$star[percentData$adjust < .001] <- "***"
    
  return (cData)
# The dataFrame of input should conform to the type below
# sample    location  origin  gene  Cp.x    location  Cp meanCP   stdDevCP ddCt    cd133
#020_N ATP5G3   B2    020_N ATP5G3 25.68         B3 25.11 25.395 0.40305087    1 negative
#020_N B2M      B4    020_N    B2M 21.34         B5 21.52 21.430 0.12727922    1 negative
#020_N CREB1    B6    020_N  CREB1 26.47         B7 26.22 26.345 0.17677670    1 negative
}

# # Run this at the end to intialise the package
# package.skeleton(name = 'qPCRcustomFunctions', 
#                  list=c('buildDataFrameForddCT', 'ddCTcalculate','extractReplicates',
#                                                       'plot_ddCt', 'splitSampleName', 'transposeLinear', 'cp', 'map'),
#                 path='~/Documents/Rscripts/', force=F) 
# #               path='/Library/Frameworks/R.framework/Versions/3.0/Resources/library/', force=F)