pkgname <- "qPCRcustomFunctions"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('qPCRcustomFunctions')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("buildDataFrameFromddCT")
### * buildDataFrameFromddCT

flush(stderr()); flush(stdout())

### Name: buildDataFrameFromddCT
### Title: Make a dataFrame that contains Cp calls and annotated sample
###   info.
### Aliases: buildDataFrameFromddCT

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
plate = transposeLinear('384wellMap.txt')
plateMap = splitSampleName(plate)
data = buildDataFrameFromddCT(plateMap, CtData)




cleanEx()
nameEx("build_ddCTmatrix")
### * build_ddCTmatrix

flush(stderr()); flush(stdout())

### Name: build_ddCTmatrix
### Title: A function to coerce ddCT values and genes into a double matrix
###   for statistical analysis
### Aliases: build_ddCTmatrix

### ** Examples

ddCtData = build_ddCTmatrix('140211_ddCtValuesCd133negPos.txt')



cleanEx()
nameEx("cp")
### * cp

flush(stderr()); flush(stdout())

### Name: cp
### Title: A example of data to manipulate
### Aliases: cp
### Keywords: datasets

### ** Examples

data(cp)
## maybe str(cp) ; plot(cp) ...



cleanEx()
nameEx("ddCTcalculate")
### * ddCTcalculate

flush(stderr()); flush(stdout())

### Name: ddCTcalculate
### Title: Calculate ddCt from raw Cp calls
### Aliases: ddCTcalculate

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
plate = transposeLinear('384wellMap.txt')
plateMap = splitSampleName(plate)
data = buildDataFrameFromddCT(plateMap, CtData)
entirePlate = c(1:384)
replicates = extractReplicates(entirePlate, data)
rawData = replicates[[3]]

rawData$ddCt = ddCTcalculate(rawData$gene, sampleOfInterest=rawData$origin,
                             houseKeepingGene='GAPDH', referenceSample='wild_type_mouse', data=rawData)
                             



cleanEx()
nameEx("extractReplicates")
### * extractReplicates

flush(stderr()); flush(stdout())

### Name: extractReplicates
### Title: Subset each Cp into its replicates
### Aliases: extractReplicates

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
plate = transposeLinear('384wellMap.txt')
plateMap = splitSampleName(plate)
data = buildDataFrameFromddCT(plateMap, CtData)

entirePlate = c(1:384)
replicates = extractReplicates(entirePlate, data)
rawData = replicates[[3]]




cleanEx()
nameEx("map")
### * map

flush(stderr()); flush(stdout())

### Name: map
### Title: An example of an experimental plate map
### Aliases: map
### Keywords: datasets

### ** Examples

data(map)
## maybe str(map) ; plot(map) ...



cleanEx()
nameEx("niceErrorBarPlot")
### * niceErrorBarPlot

flush(stderr()); flush(stdout())

### Name: niceErrorBarPlot
### Title: Generates a nice looking barplot that has error bars
### Aliases: niceErrorBarPlot

### ** Examples

There is no working example yet. Need to work on this.



cleanEx()
nameEx("niceGroupedBarPlot")
### * niceGroupedBarPlot

flush(stderr()); flush(stdout())

### Name: niceGroupedBarPlot
### Title: Make a nice barchart to summarise a qPCR experiment.
### Aliases: niceGroupedBarPlot

### ** Examples

This function needs more work to be of use. I think it is parsing of the column names



cleanEx()
nameEx("qPCRcustomFunctions-package")
### * qPCRcustomFunctions-package

flush(stderr()); flush(stdout())

### Name: qPCRcustomFunctions-package
### Title: This package contains conveinace functions to analyse qPCR
###   generated using 384 well plates, specifically the lightcycler 480
### Aliases: qPCRcustomFunctions-package qPCRcustomFunctions
### Keywords: package

### ** Examples

~~ simple examples of the most important functions ~~
plate = transposeLinear('384wellMap.txt')
plateMap = splitSampleName(plate)
data = buildDataFrameFromddCT(plateMap, CtData)

entirePlate = c(1:384)
replicates = extractReplicates(entirePlate, data)
rawData = replicates[[3]]

rawData$ddCt = ddCTcalculate(rawData$gene, sampleOfInterest=rawData$origin,
                             houseKeepingGene='GAPDH', referenceSample='wild_type_mouse', data=rawData)
                             
p = plot_ddCt(ddCt)~origin, data, 'Some data to plot', yaxisLabel='ddCT')                         



cleanEx()
nameEx("splitSampleName")
### * splitSampleName

flush(stderr()); flush(stdout())

### Name: splitSampleName
### Title: Split a column containing origin and gene into 2 columns
### Aliases: splitSampleName

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
plate = transposeLinear('textFile384well.txt')
data = splitSampleName(plateMap) 



cleanEx()
nameEx("summariseStatistics_ddCt")
### * summariseStatistics_ddCt

flush(stderr()); flush(stdout())

### Name: summariseStatistics_ddCt
### Title: Generate mean and standard error.
### Aliases: summariseStatistics_ddCt

### ** Examples

# Bind individual biological replicates into 1 dataframe
cd133negPos = rbind(cd133_20, cd133_30a, cd133_41)

# Populate a column listing whether the sample is CD133+ or CD133-. The summariseStatistics_ddCt function uses this
cd133negPos$cd133 = ifelse(grepl('*_N', cd133negPos$origin.x), 'negative', 'positive')

# Return summary statistics
cData = summariseStatistics_ddCt(cd133negPos)



cleanEx()
nameEx("transposeLinear")
### * transposeLinear

flush(stderr()); flush(stdout())

### Name: transposeLinear
### Title: Convert plate map to R object
### Aliases: transposeLinear

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
read.delim('well384Map.txt')
function (well384Map, linearMapFile = "output.txt") 



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
