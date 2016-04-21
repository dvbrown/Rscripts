build_ddCTmatrix <-
function(ddCtFile, originColumn=2, geneColumn=3, ddCtColumn=9, output='matrix.txt') {
    # A function to coerce ddCT values and genes into a double matrix for statistical analysis
    # origin is the name of the sample eg #020 CD133 negative
    pythonCall = paste('../exec/buildNumericMatrix_ddCt.py', '-i', ddCtFile,
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
