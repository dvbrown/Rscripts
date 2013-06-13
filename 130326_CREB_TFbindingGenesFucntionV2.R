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
genes = read.table('130326_promoterGeneTester.txt',header=T, stringsAsFactors=F)
genes = genes$Symbol

#compare the 2 motifs for the transcription factor
pfm.creb1.hPDI = new('pfm', mat=query(MotifDb, transcriptionFactor)[[1]], name="database 1")
pfm.creb1.jaspar <- new("pfm", mat=query(MotifDb, transcriptionFactor)[[2]], name="database 2")
plotMotifLogoStack(DNAmotifAlignment(c(pfm.creb1.hPDI, pfm.creb1.jaspar)))

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

countMotifInPromoter <- function (genePromoter, posCountMatrix) {
  #First argument is the promoter you want to scan. Second is the postion count matrix of the transciptionFactor
  #TO DO enclose in a try except block in case it can't find the gene name
  orfs = getOrf(genes)
  #Store the genomic coordinates in a GR object and extract the promter sequences
  grl = transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, vals=list(gene_id=orfs))
  promoter.seqs = getPromoterSeq(grl, BSgenome.Hsapiens.UCSC.hg19, upstream=1000, downstream=0)
  #reduce the list structure of the DNA string objects
  promoter.seqs = unlist(promoter.seqs[1])
  
  #now the matching. The indexing of the promoter list means match only the first transcript ID
  pwm.hits = matchPWM(pwm=posCountMatrix, subject=promoter.seqs, min.score="90%")
  print(pwm.hits)
  count = countPWM(pwm=posCountMatrix, subject=promoter.seqs, min.score="90%")
  return (count)
}

getOrf = function(genes) {
  orfs = as.character(mget(genes, org.Hs.egALIAS2EG))
  #if(class(orfs) == "try-error") return(count = 0);
  return (orfs)
}

pcmTranscriptionFactor = getPositionCountMatrix(transcriptionFactor, 2)

for (gene in genes) {
  print(gene)
  promtoterMotifCounts = try(countMotifInPromoter(genePromoter=gene,posCountMatrix=pcmTranscriptionFactor));
  if(class(promtoterMotifCounts) == 'try-error') next;
}


counts = cbind(genes[,], promtoterMotifCounts)
counts
write.table(counts, './130326_motifOut.txt', sep='\t')
#return to old directory
setwd(currentDir)