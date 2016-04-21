library(GenomicRanges)
data(cpgIslands) #load example data
class(cpgIslands)

chr <- as.character(unique(seqnames(cpgIslands)))
gen <- genome(cpgIslands) #retreive sequence coordinates
atrack <- AnnotationTrack(cpgIslands, name = "CpG") #build the track object
gtrack <- GenomeAxisTrack()
#Since a GenomeAxisTrack object is always relative to the other tracks that are plotted, there is little need for 
#additional arguments. Essentially, the object just tells the plotTracks function to add a genomic axis to the plot. 
#Nonetheless, it represent a separate annotation track just as the CpG island track does. 
#We can pass this additional track on to plotTracks in the form of a list.
plotTracks(list(gtrack, atrack))
itrack <- IdeogramTrack(genome = gen, chromosome = chr) #Get the chromosome image from UCSC
#Similar to the previous examples, we stick the additional track object into a list in order to plot it.
plotTracks(list(itrack, gtrack, atrack)) 

#build a gene model track
data(geneModels)
grtrack <- GeneRegionTrack(geneModels, genome = gen,
                             chromosome = chr, name = "Gene Model")
plotTracks(list(itrack, gtrack, atrack, grtrack))
#zoom in on a particular reigon using the from and to arguments
plotTracks(list(itrack, gtrack, atrack, grtrack), from = 2.5e+07, to = 2.8e+07)

library(BSgenome.Hsapiens.UCSC.hg19)
strack <- SequenceTrack(Hsapiens, chromosome = chr) #get the actual sequence
plotTracks(list(itrack, gtrack, atrack, grtrack, strack), from = 26450430, to = 26450490, cex = 0.8)

set.seed(255)
lim <- c(26463500, 26495000)
coords <- sort(c(lim[1], sample(seq(from = lim[1], to = lim[2]), 99), lim[2]))
dat <- runif(100, min = -10, max = 10) #make some dummy data
#The individual rows in a numeric matrix are considered to be different data groups or samples, and the 
#columns are the raster intervals in the genomic coordinates.
dtrack <- DataTrack(data = dat, start = coords[-length(coords)], end = coords[-1], chromosome = chr, genome = gen, name = "Uniform")
plotTracks(list(itrack, gtrack, atrack, grtrack, dtrack), from = lim[1], to = lim[2])
#change the dot plot to a histogram in the data object
plotTracks(list(itrack, gtrack, atrack, grtrack, dtrack), from = lim[1], to = lim[2], type = "histogram")

#example usage of additional arguments
grtrack <- GeneRegionTrack(geneModels, genome = gen, chromosome = chr, name = "Gene Model", showId = TRUE,
                           background.title = "brown")
displayPars(grtrack) <- list(background.panel = "#FFFEDB")
plotTracks(list(itrack, gtrack, atrack, grtrack),
           from = lim[1], to = lim[2])