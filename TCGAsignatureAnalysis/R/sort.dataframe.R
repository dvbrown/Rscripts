sort.dataframe <-
function(dataframe, columnNumber, highFirst=TRUE) {
  #if highFirst is true the return the list from high to low
  new.dataframe = dataframe[with(dataframe, order(dataframe[,columnNumber], decreasing=highFirst)),]
  return (new.dataframe)
}
