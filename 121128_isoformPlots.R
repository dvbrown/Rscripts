#try to read in some gapped alignments. Uses the GViz library
library(ggbio) 
library(IRanges)
library(GenomicRanges)
library(Rsamtools)
library(GenomicFeatures)
library(cummeRbund)
setwd('~/Documents/RNAdata/RNAseqAnalysis/121116_cuffOutFull/')
#Intialise the connection to the SQL lite database. Do not rebuild the database.
cuff = readCufflinks(dir=getwd(), gtfFile=('~/Documents/public-datasets/annotationFiles/merged.gtf'), genome='hg19')
txdb=makeTranscriptDbFromUCSC(genome='hg19',tablename='ensGene')
#extract the gene of interest by modifiying gene
myGene = getGene(cuff, 'gene')
myChr<-unique(features(myGene)$seqnames)
myStart<-min(features(myGene)$start)
myEnd<-max(features(myGene)$end)
genome<-'hg19'

extractOverlaps = function(genomicRanges) { #take a genomic ranges object that was subset for a particular gene. Doesn't work as function.
  #extract the exon specific coordinates from ENsembl
  tx_by_exon=transcriptsBy(txdb,by=c('exon', 'cds'))
  tx = subsetByOverlaps(tx_by_exon, genomicRanges, type='any')
  result = countOverlaps(genomicRanges, tx, type='any') #check this is the right way around
  #return data for constructing the data track in GViz
  return (result)
}

getGeneAlignments = function(pathToBam) { #take a string represing a bam file and return 
  gene_coords <- GRanges(seqnames = myChr, ranges = IRanges(start = myStart, end = myEnd))
  bam_fields <- c("pos", "qwidth")
  param <- ScanBamParam(which = gene_coords, what = bam_fields)
  bam_file <- pathToBam
  bam <- readGappedAlignments(bam_file, param = param, use.names=T)
  gr <- granges(bam)
  seqlevels(gr)=gsub("^","chr",seqlevels(gr)) #match up the chromsome names from Bam file to Ensembl data
  return (gr)
}

#Plot a gene if the gene region has any coverage
gene_plot <- if(length(width(gr)) >= 1){
  #exons contains the start and stop regions for each exon in the gene
  exons <- data.frame(xmin = myStart, xmax = myEnd, ymin = 0, ymax = Inf)
  cov_plot <- autoplot (
    gr, 
    stat = "coverage",
    geom = "area",
    colour = "black",
    fill = "green",
    xlab = "Position",
    ylab = "Coverage",
    main = "Coverage across gene" 
    #Marks the exon regions out so that you can see how much of the sequencing is on target  			
  ) + geom_rect(
    data = exons, 
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
    fill = "yellow", 
    colour = "black", alpha = 0.5
  )
}
#plot the graph				
gene_plot

#make a trackList and build it up with tracks
trackList<-list()
ideoTrack <- IdeogramTrack(genome = genome, chromosome = myChr) #loads the chromosome where myGene resides
genetrack<-makeGeneRegionTrack(myGene) #loads the gene model for myGene
dTrack = DataTrack(range=gr, genome='hg19', data=x, type='hist', window=200, name='aligned reads') # the window argument is the number of bins
trackList = c(ideoTrack, genetrack, dTrack)
plotTracks(trackList,from=myStart-1000,to=myEnd+1000)