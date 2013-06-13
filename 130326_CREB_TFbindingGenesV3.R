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

transcriptionFactor = as.character('CREB1')
genes = as.character('LOC90410')

#search the database for the transcription factor of choice
query(MotifDb, transcriptionFactor)
#the representation of the TF binding sequence
pfm.creb1.hPDI = query(MotifDb, transcriptionFactor )[[1]]
pfm.creb1.jaspar = query(MotifDb, transcriptionFactor)[[2]]
#turn the pfm into a count matrix by multiplying by 100
pcm.creb1.hPDI = round(100 * pfm.creb1.hPDI)
pcm.creb1.jaspar = round(100 * pfm.creb1.jaspar)


#compare the 2 motifs for the transcription factor
pfm.creb1.hPDI = new('pfm', mat=query(MotifDb, transcriptionFactor)[[1]], name="CREB1-hPDI")
pfm.creb1.jaspar <- new("pfm", mat=query(MotifDb, transcriptionFactor)[[2]], name="CREB-jaspar")
plotMotifLogoStack(DNAmotifAlignment(c(pfm.creb1.hPDI, pfm.creb1.jaspar)))



#retreive the ORF IDs.
#TO DO enclose in a try except block in case it can't find the gene name
getOrf = function(genes) {
  orfs = try(as.character(mget(genes, org.Hs.egALIAS2EG)));
  if(class(orfs) == "try-error") next;
return (orfs)
}
orfs = getOrf(genes)
#orfs = as.character(mget(genes, org.Hs.egALIAS2EG))
#Store the genomic coordinates in a GR object and extract the promter sequences
grl = transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, vals=list(gene_id=orfs))
promoter.seqs = getPromoterSeq(grl, BSgenome.Hsapiens.UCSC.hg19, upstream=1000, downstream=0)

#reduce the list structure of the DNA string objects
promoter.seqs = unlist(promoter.seqs[1])

#now the matching. The indexing of the promoter list means match only the first transcript ID
pwm.hits = matchPWM(pwm=pcm.creb1.jaspar, subject=promoter.seqs, min.score="90%")
count = countPWM(pwm=pcm.creb1.jaspar, subject=promoter.seqs, min.score="90%")
count
pwm.hits
#A function to append matches into a vector

#return to old directory
setwd(currentDir)