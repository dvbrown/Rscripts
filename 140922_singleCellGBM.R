setwd('~/Documents/public-datasets/RNA-seq/anoop2014_singleCellGBM/')

subsetSamples <- function (dataFrame, sampleNameStub) {
    # The dataframe with the measurements
    # A string representing the basename of the sample eg MGH26
    sample = paste(sampleNameStub, '*', sep='_')
    samples = grep(sample, colnames(dataFrame), value=T)
    result = dataFrame[,samples]
    return (result)
}


data = read.delim('GSE57872_GBM_data_matrix.txt', row.names=1)
annotation = read.delim('sample.txt')

# Extract the interesting samples
interesting = c('MGH26', 'MGH28', 'MGH29', 'MGH30', 'MGH31')
mgh26Data = subsetSamples(data, 'MGH26')
mgh28Data = subsetSamples(data, 'MGH28')
mgh29Data = subsetSamples(data, 'MGH29')
mgh30Data = subsetSamples(data, 'MGH30')
mgh31Data = subsetSamples(data, 'MGH31')
