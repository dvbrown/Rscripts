#!/usr/bin/env Rscript

library(optparse)

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")
option_list <- list(
    make_option(c("-a", "--about"), action="store_true", default=FALSE,
                help="An Rscript to find candidate binding sites for known transcription factors via seq matching (Bioconductor). Based on the tutorial http://www.bioconductor.org/help/workflows/gene-regulation-tfbs/"),
    make_option(c("-t", "--transcriptionFactor"), action="store", type='character', default='CREB1',
                help="The gene symbol of the transcription factor eg CREB1"),
    make_option(c("-g", "--geneList"), action="store",type = 'character', default='~/Documents/CREB/transcriptionFactorAnalysis/130325_TFCREB/130326_promoterGeneTester.txt',
                help="A tab delimited text file (genes as column 1) of listing the geneSymbols of the promoters you wish to interrogate"),
    make_option(c("-o", "--outFile"), action="store", type='character', default='output.txt',
                help="The file you wish to output results as a tab delimited text file"),
    make_option(c("-u", "--upstream"), action="store", type='integer', default=1000,
                help="The distance upstream of the transcription start site you wish to scan"),
    make_option(c("-s", "--scanDownstream"), action="store", type='integer', default=0,
                help="The distance downstream of the transcription start site you wish to scan")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))
###########################################################################################################
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
source('/Users/d.brown6/Documents/Rscripts/transcriptionFactorMotif/130326_transcriptionFactorMatchingFunc.R')

currentDir = getwd()
setwd('/Users/d.brown6/Documents/CREB/transcriptionFactorAnalysis/130325_TFCREB/')

#global variables parsed from command line
transcriptionFactor = opt$transcriptionFactor
geneFile = opt$geneList
outFile = opt$outFile
genes = as.matrix(read.table(geneFile, header=T, stringsAsFactors=F))
genes = unique(genes)
dataBaseToUse = opt$databaseNo
upStream = opt$upstream
downStream = opt$scanDownstream

#The function calls
orfs = geneNameToID(genes)
pcmCREBfull = getPositionCountMatrix(transcriptionFactor, 2)
pcmCREBhalf = getPositionCountMatrix(transcriptionFactor, 1)

#Use apply to use the countMotifInPromoter function across the vector of geneIDs
CREBfullMotifCounts = apply(orfs, 1, countMotifInPromoter, posCountMatrix=pcmCREBfull, up=upStream,down=downStream)
CREBhalfMotifCounts = apply(orfs, 1, countMotifInPromoter, posCountMatrix=pcmCREBhalf, up=upStream,down=downStream)

counts = cbind(orfs[,], CREBfullMotifCounts, CREBhalfMotifCounts)
colnames(counts) = c('enterezID', 'fullCREBsites', 'halfCREBsites')

#print(counts)
write.table(counts, outFile, sep='\t')

#return to old directory
setwd(currentDir)