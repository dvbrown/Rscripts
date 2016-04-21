#this is currently the best method
library(ggbio)
library(GenomicRanges)
library(GenomicFeatures)
library(Rsamtools)
library(cummeRbund)
setwd('~/Documents/RNAdata/RNAseqAnalysis/121116_cuffOutFull/')
library(BSgenome.Hsapiens.UCSC.hg19)
#build the gene model plots for your gene of choice
txdb=makeTranscriptDbFromUCSC(genome='hg19',tablename='ensGene')
data(genesymbol, package = 'biovizBase')

gene = "ANXA1" #change this gene depending on what you want to see

geneModelFull = autoplot(txdb, which = genesymbol[gene], names.expr = "tx_name")
geneModelReduce = autoplot(txdb, which = genesymbol[gene], stat = "reduce", color = "brown", fill = "brown")

#Intialise the connection to the SQL lite database. Do not rebuild the database. Use these figures to scan the bam file
cuff = readCufflinks(dir=getwd(), gtfFile=('~/Documents/public-datasets/annotationFiles/merged.gtf'), genome='hg19')
myGene = getGene(cuff, gene)
myChr<-unique(cummeRbund::features(myGene)$seqnames)
myStart<-min(cummeRbund::features(myGene)$start)
myEnd<-max(cummeRbund::features(myGene)$end)
genome<-'hg19'

#read in bam file
gene_coords <- GRanges(seqnames = myChr, ranges = IRanges(start = myStart, end = myEnd))
bam_fields <- c("pos", "qwidth")
param <- ScanBamParam(which = gene_coords, what = bam_fields)
bam_file1 <- "~/Documents/RNAdata/RNAseqAnalysis/121121_TopHatNovelAlignment/CD133nPoolSorted.bam" #put the bam file in here
bam_file2 <- "~/Documents/RNAdata/RNAseqAnalysis/121121_TopHatNovelAlignment/CD133pPoolSorted.bam"
bam1 <- readGappedAlignments(bam_file1, param = param, use.names=T) #only read the bam file where your gene is
bam2 <- readGappedAlignments(bam_file2, param = param, use.names=T)
grn <- granges(bam1) #this is the bam file subset by gene name
grp <- granges(bam2)
#add sample annotation to Granges object
values(grn)$sample = 'negative'
values(grp)$sample = 'positive'
#concatentate the Granges objects
comparison = c(grn, grp)
myChr=gsub("^","chr",myChr)
ideogram = plotIdeogram(genome='hg19', subchr=myChr) #must use chr notation. Look at value of myChr and append this here
compPlot = autoplot(comparison, geom='line', stat='coverage', facets=sample~seqnames)
tracks(chr=ideogram, complete=geneModelFull, summary=geneModelReduce, data=compPlot, heights=c(2,3,2,5), main= 'ANXA1 gene coverage')

geneCoverage = autoplot(grp, geom='area', stat='coverage')
#Plot the objects together
tracks(chr=ideogram, complete=geneModelFull, summary=geneModelReduce, data=geneCoverage, heights=c(1,3,2,5), main='TMSB4X gene coverage')
linear = plotGrandLinear(gr, coord='genome', geom='point', stat='coverage', space.skip=0.01, color=gene_coords) 
plotRangesLinkedToData(gr) #these plots need the whole bam file.