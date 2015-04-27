bindSignatureToSurvival <-
function(signatureScore, clinical) {
    #take the annotated clinical data and bind it a computed gene signature score for analysis    
    clinical = clinical[,c(1,4)] #get just the survival and censorship. FIX DIS
    geneSet = as.data.frame(signatureScore)
    set = intersect(row.names(geneSet), row.names(clinical)) #get the overlap
    data = cbind(geneSet[set,], clinical[set,]) #bind the 2 datasets together based on overlap
    colnames(data) = c('sigScore', 'survival', 'censorship')
    return(data)
}
