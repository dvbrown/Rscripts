#Hold the functions to be imported into main

geneNameToID <- function (genes) {
  #Takes a matrix of gene Symbols as argument and returns the geneID mappings.
  holder = vector(mode = 'numeric', len=length(genes))
  i = 1
  #check that gene names match a gene ID
  for (gene in genes) {
    print(gene)
    x = try(as.character(mget(gene, org.Hs.egALIAS2EG)))
    if(class(x) == "try-error") {next} 
    else {holder[i] = x}
    i = i + 1
  }
  #retain only those geneIDs that have a match (ie are not 0)
  orfs = as.matrix(holder[which(holder!='0')])
  return (orfs)
}

#search the database for the transcription factor of choice

getPositionCountMatrix <- function (transcriptionFactor, databaseToUse) {
  #First argument is a string of the transcription factor. 
  #Second is the integer of database number you wish to use. ie '2' in the case of CREB
  query(MotifDb, transcriptionFactor)
  #the representation of the TF binding sequence
  pfmTransfactor = query(MotifDb, transcriptionFactor )[[databaseToUse]]
  #turn the pfm into a count matrix by multiplying by 100
  pcmTransFactor = round(100 * pfmTransfactor)
  return (pcmTransFactor)
}

countMotifInPromoter <- function (orfs, posCountMatrix, up, down) {
  #First argument is the promoter you want to scan. Second is the postion count matrix of the transciptionFactor
  #internal function call to get geneID from gene symbols
  #Store the genomic coordinates in a GR object and extract the promter sequences
  print(orfs)
  grl = transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, vals=list(gene_id=orfs))
  promoter.seqs = getPromoterSeq(grl, BSgenome.Hsapiens.UCSC.hg19, upstream=up, downstream=down)
  #reduce the list structure of the DNA string objects
  #promoter.seqs = unlist(promoter.seqs[1])
  
  promoter.seqs = try(unlist(promoter.seqs[1]))
  if (class(promoter.seqs) == "try-error") 
  {count = NA}
  #now the matching. The indexing of the promoter list means match only the first transcript ID
  else {pwm.hits = matchPWM(pwm=posCountMatrix, subject=promoter.seqs, min.score="90%")
        #print(pwm.hits)
        count = countPWM(pwm=posCountMatrix, subject=promoter.seqs, min.score="90%")}
  return (count)
}

transcriptIDtoSymbol <- function (bedFile, key) {
    #Takes as input the bed file output by the nearest BEDtools program in BED tools
    transcriptID = data[,13]
    mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    
    x = getBM(
        filters= bioMartKey, 
        attributes= c("ensembl_transcript_id","ensembl_gene_id", "external_gene_id", "entrezgene", "description"),
        values= transcriptID,
        mart= mart)
    result = merge(data, x, by.x='V13', by.y='ensembl_transcript_id')
    result = result[,c(1,2,3,4,5,6,7,10,11,12,13,22,23,24,25)]
    return (result)
}

#This next function needs to be fixed. It is not yet called in main
enterzIDtoGeneSymbol <- function () {
  #bind the vector of motif matches to the geneSymbol
  geneSymbol = orfs[,]
  x = org.Hs.egSYMBOL
  mappedGenes = mappedkeys(x)
  xx = as.list(x[mappedGenes])
  symbol = xx[geneSymbol]
  for (i in symbol) {
      symbol = unlist(symbol)
      return (symbol)
    }
}