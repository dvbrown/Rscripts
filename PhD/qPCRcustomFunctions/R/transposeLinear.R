transposeLinear <-
function(well384Map, linearMapFile='output.txt') {
  pythonCall = paste('../exec/transposeLinear.py -i', well384Map, '>', linearMapFile, sep=' ')
  system(pythonCall)
  sampleLabels = read.delim(linearMapFile, header=F)
  colnames(sampleLabels) = c('location', 'sample')
  return (sampleLabels)
}
