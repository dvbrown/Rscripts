summariseStatistics_ddCt <-
function (dataFrame, groupVariableA='cd133', groupVariableB='gene.x') {
    # Generate N, mean, sd and se statistics for a dataframe
    cData = ddply(cd133negPos, c(groupVariableA, groupVariableB), summarise,
                N    = sum(!is.na(ddCt)),
                mean = mean(ddCt, na.rm=TRUE),
                sd   = sd(ddCt, na.rm=TRUE),
                se   = sd / sqrt(N) )
  return (cData)
# The dataFrame of input should conform to the type below
# sample    location  origin  gene  Cp.x    location  Cp meanCP   stdDevCP ddCt    cd133
#020_N ATP5G3   B2    020_N ATP5G3 25.68         B3 25.11 25.395 0.40305087    1 negative
#020_N B2M      B4    020_N    B2M 21.34         B5 21.52 21.430 0.12727922    1 negative
#020_N CREB1    B6    020_N  CREB1 26.47         B7 26.22 26.345 0.17677670    1 negative
}
