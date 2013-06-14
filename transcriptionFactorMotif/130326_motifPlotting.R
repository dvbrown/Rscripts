#!/usr/bin/env Rscript

library(optparse)
library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(org.Hs.eg.db)

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")
option_list <- list(
    make_option(c("-h", "--help"), action="store_true", default=FALSE,
                help="An Rscript to plot the transcription factor motifs to find out which to use")
    make_option(c("-t", "--transcriptionFactor"), action="store_true", default=FALSE,
                help="The gene symbol of the transcription factor eg CREB1")
)

transcriptionFactor = as.character('CREB1')

#compare the 2 motifs for the transcription factor
pfm.creb1.hPDI = new('pfm', mat=query(MotifDb, transcriptionFactor)[[1]], name="Human Protein-DNA Interactome")
pfm.creb1.jaspar <- new("pfm", mat=query(MotifDb, transcriptionFactor)[[2]], name="Jaspar")
plotMotifLogoStack(DNAmotifAlignment(c(pfm.creb1.hPDI, pfm.creb1.jaspar)))