#!/usr/bin/env Rscript
#An Rscript to find candidate binding sites for known transcription factors via seq matching (Bioconductor)
#Based on the tutorial http://www.bioconductor.org/help/workflows/gene-regulation-tfbs/
#To match motifs in a promoter, these steps are required:
#Retrieve the binding motif (the position frequency matrix, or PFM) of a given transcription factor
#Retrieve the promoter regions for a set of candidate targets
#Identify the sequence matches of the binding motif in the the genes' promoter regions

currentDir = getwd()
setwd('/Users/d.brown6/Documents/CREB/transcriptionFactorAnalysis/130325_TFCREB/')

transcriptionFactor = 'CREB1'
genes = as.character(c('GPR3','SPAG4', 'STX10', 'ACBD4'))

library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

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
orfs = as.character(mget(genes, org.Hs.egALIAS2EG))
#Store the genomic coordinates in a GR object and extract the promter sequences
grl = transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')[orfs]
promoter.seqs = getPromoterSeq(grl, BSgenome.Hsapiens.UCSC.hg19, upstream=1000, downstream=0)

#reduce the list structure of the DNA string objects
print(class(promoter.seqs))
promoter.seqs = unlist(promoter.seqs)
print(class(promoter.seqs))

#now the matching. This is not working. encapsulate the promoter finding per gene as a function.
for (sequence in promoter.seqs) {
    pwm.hits = sapply(subject=sequence, function(pseq)matchPWM(pwm=pcm.creb1.jaspar, pseq, min.score="90%"))
    fos.jaspar = length(pwm.hits) 
    }
fos.jaspar

#return to old directory
setwd(currentDir)