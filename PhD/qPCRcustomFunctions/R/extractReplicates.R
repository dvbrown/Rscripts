extractReplicates <-
function (indexes, ctData) {
  # Retreive the indexes of the 384 wellplate
  indexes = c(1:384)
  
  # Keep only cases with data in them as the merge function doesn't work with NAs
  CtData = data[complete.cases(data[,3]),]
  # Subset each Cp into its replicates. Takes a vector with the indexes to to subset and then takes the
  # even entries and odd entries separately from the dataframe containing cp values
  even = indexes[indexes%%2 == 0]
  odd = indexes[indexes%%2 == 1]
  
  rep1 = CtData[odd, c(1:6)]
  rep1 = na.omit(rep1)
  rep2 = CtData[even, c(1:6)]
  rep2 = na.omit(rep2)
  boundData = merge(rep1, rep2, by.x='sample', by.y='sample')
  ################ Remove columns that do not add information
  usefulData = boundData[,c(1,2,3,4,6,7,11)]
  # Compute the mean and the standard deviation of the replicates
  usefulData$meanCP = rowMeans(cbind(usefulData$Cp.x, usefulData$Cp.y))
  usefulData$stdDevCP = apply(cbind(usefulData$Cp.x, usefulData$Cp.y), 1, sd)
  # Package the output in a list
  result = list(rep1, rep2, usefulData)
  return (result)
}
