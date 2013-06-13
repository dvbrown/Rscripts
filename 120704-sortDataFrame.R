sort.dataframe <- function(dataframe, column, highFirst=TRUE) {
  #if highFirst is true the return the list from high to low
  new.dataframe = dataframe[with(dataframe, order(column, decreasing=highFirst)),]
  new.dataframe
}