#!/usr/bin/env Rscript
###########################################################################################################
#An Rscript to find candidate binding sites for known transcription factors via seq matching (Bioconductor)
#Based on the tutorial http://www.bioconductor.org/help/workflows/gene-regulation-tfbs/
#To match motifs in a promoter, these steps are required:
#Retrieve the binding motif (the position frequency matrix, or PFM) of a given transcription factor
#Retrieve the promoter regions for a set of candidate targets
#Identify the sequence matches of the binding motif in the the genes' promoter regions
###########################################################################################################

library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
currentDir = getwd()
setwd('/Users/d.brown6/Documents/CREB/transcriptionFactorAnalysis/130325_TFCREB/')

#global variables parsed from command line
transcriptionFactor = as.character('CREB1')
genes = as.matrix(read.table('130326_promoterGeneTester.txt',header=T, stringsAsFactors=F))

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

countMotifInPromoter <- function (orfs, posCountMatrix) {
  #First argument is the promoter you want to scan. Second is the postion count matrix of the transciptionFactor
  #internal function call to get geneID from gene symbols
  #Store the genomic coordinates in a GR object and extract the promter sequences
  print(orfs)
  grl = transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, vals=list(gene_id=orfs))
  promoter.seqs = getPromoterSeq(grl, BSgenome.Hsapiens.UCSC.hg19, upstream=1000, downstream=0)
  #reduce the list structure of the DNA string objects
  #promoter.seqs = unlist(promoter.seqs[1])
  
  promoter.seqs = try(unlist(promoter.seqs[1]))
  if (class(promoter.seqs) == "try-error") 
    {count = NA}
  #now the matching. The indexing of the promoter list means match only the first transcript ID
  else {pwm.hits = matchPWM(pwm=posCountMatrix, subject=promoter.seqs, min.score="90%")
    print(pwm.hits)
    count = countPWM(pwm=posCountMatrix, subject=promoter.seqs, min.score="90%")}
  return (count)
}

#compare the 2 motifs for the transcription factor
#pfm.creb1.hPDI = new('pfm', mat=query(MotifDb, transcriptionFactor)[[1]], name="database 1")
#pfm.creb1.jaspar <- new("pfm", mat=query(MotifDb, transcriptionFactor)[[2]], name="database 2")
#plotMotifLogoStack(DNAmotifAlignment(c(pfm.creb1.hPDI, pfm.creb1.jaspar)

#The function calls
orfs = geneNameToID(genes)
pcmTranscriptionFactor = getPositionCountMatrix(transcriptionFactor, 2)

#Use apply to use the countMotifInPromoter function across the vector of geneIDs
promtoterMotifCounts = apply(orfs, 1, countMotifInPromoter, posCountMatrix=pcmTranscriptionFactor)

#bind the vector of motif matches to the geneSymbol
counts = cbind(orfs[,], promtoterMotifCounts)
colnames(counts) = c('enterezID', 'motifMatches')
print(counts)
write.table(counts, './130326_motifOut.txt', sep='\t')

#return to old directory
setwd(currentDir)