library(cummeRbund)
setwd('~/Documents/RNAdata/RNAseqAnalysis/121116_cuffOutFull/')
#Intialise the connection to the SQL lite database. Do not rebuild the database.
cuff = readCufflinks(dir=getwd(), gtfFile=('~/Documents/public-datasets/annotationFiles/merged.gtf'), genome='hg19')
disp = dispersionPlot(genes(cuff))
genes.scv<-fpkmSCVPlot(genes(cuff))
isoforms.scv<-fpkmSCVPlot(isoforms(cuff))
brep<-csBoxplot(genes(cuff),replicates=T)
s<-csScatterMatrix(genes(cuff))
dend.rep<-csDendro(genes(cuff),replicates=T)

myGene = getGene(cuff, 'EGF')
head(features(myGene))
trackList<-list()
myStart<-min(features(myGene)$start)
myEnd<-max(features(myGene)$end)
myChr<-unique(features(myGene)$seqnames)
genome<-'hg19'
ideoTrack <- IdeogramTrack(genome = genome, chromosome = myChr) #loads the chromosome where myGene resides
trackList<-c(trackList,ideoTrack)
axtrack<-GenomeAxisTrack()
trackList<-c(trackList,axtrack)
genetrack<-makeGeneRegionTrack(myGene) #loads the gene model for myGene
trackList<-c(trackList,genetrack)
biomTrack<-BiomartGeneRegionTrack(genome=genome,chromosome=as.character(myChr), + start=myStart,end=myEnd,name="ENSEMBL",showId=T)
trackList<-c(trackList,biomTrack)
dtrack <- DataTrack(data = myGene, start = coords[-length(coords)], end = coords[-1], 
                    chromosome = chr, genome = genome, name = "Uniform")
trackList<-c(trackList,conservation)
plotTracks(trackList,from=myStart-2000,to=myEnd+2000)

#try to read in some gapped alignments
library(ggbio) 
library(IRanges)
library(GenomicRanges)
library(Rsamtools)
library(GenomicFeatures)

countreads = function(GappedAlignments) { #Take a gapped alignments object (Bam file) and return a vector of counts.
  #change the bam file chromosome names to match the genome annotation, ie the tx_by_gene object
  seqlevels(gr)=gsub("^","chr",seqlevels(gr))
  reads=GRanges(seqnames=new_read_chr_names,ranges=IRanges(start=myStart,end=myEnd), strand=rep("*",length(gr)))
  #extract read counts
  counts=countOverlaps(tx_by_exon,reads)
  return (counts)
}

myGene = getGene(cuff, 'EGF')
myChr<-unique(features(myGene)$seqnames)
myStart<-min(features(myGene)$start)
myEnd<-max(features(myGene)$end)
genome<-'hg19'

txdb=makeTranscriptDbFromUCSC(genome='hg19',tablename='ensGene')
#extract the exon specific coordinates from ENsembl
tx_by_exon=transcriptsBy(txdb,by-c('exon', 'cds'))
tx = subsetByOverlaps(tx_by_exon, gr)
x = countOverlaps(gr, tx)

gene_coords <- GRanges(seqnames = 4, ranges = IRanges(start = myStart, end = myEnd))
bam_fields <- c("pos", "qwidth")
param <- ScanBamParam(which = gene_coords, what = bam_fields)
bam_file <- "~/Documents/RNAdata/RNAseqAnalysis/121121_TopHatNovelAlignment/CD133nPoolSorted.bam"
bam <- readGappedAlignments(bam_file, param = param, use.names=T)
gr <- granges(bam)
seqlevels(gr)=gsub("^","chr",seqlevels(gr))

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
#make a trackList and plot it
trackList<-list()
ideoTrack <- IdeogramTrack(genome = genome, chromosome = myChr) #loads the chromosome where myGene resides
genetrack<-makeGeneRegionTrack(myGene) #loads the gene model for myGene
#build the overlaps object
dTrack = DataTrack(range=gr, genome='hg19', data=x)
trackList = c(ideoTrack, genetrack, dTrack)
plotTracks(trackList,from=myStart-2000,to=myEnd+2000)