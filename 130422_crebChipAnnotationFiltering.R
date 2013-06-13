library(BSgenome.Hsapiens.UCSC.hg19)
library(ggbio)
library(GenomicRanges)
library(GenomicFeatures)
library(hexbin)
setwd('~/Documents/CREB/ChIPseqENCODE/insectBED1000/')

data = read.delim('130422_mergedOutput.annotateOneExpt.countMotifs.txt', row.names=1, stringsAsFactors = F)
data = data[,1:21]

#convert to right type
data$PeakScore = as.integer(data$PeakScore)
data$Distance.to.TSS = as.integer(data$Distance.to.TSS)
data$intervalSize = data$End - data$Start
#normalise the tag and motif number for the size of the interval
data$tagsPerInterval = (data$PeakScore / data$intervalSize )
data$motifPerInterval = (data$CREBmotifNo / data$intervalSize)

#the filters
dataPeakScore = subset(data, data$PeakScore >= 12)
dataTSSdist = subset(dataPeakScore, dataPeakScore$Distance.to.TSS <= 3000 & dataPeakScore$Distance.to.TSS >= -500)
dataMotif = subset(dataTSSdist, dataTSSdist$CREBmotifNo > 0)
dataMotifDensity = subset(dataMotif, dataMotif$motifPerInterval <= 1000)
#dataProtein = subset(dataMotifDensity, dataMotifDensity$Gene.Type == 'protein-coding')

#make some graphs
par(mfrow=c(2,1))
plot(data$Distance.to.TSS, data$motifPerInterval, type='h', xlab='Distance to transcriptional start site',
     ylab='Motifs per base pair', main='Density of CREB motifs', xlim=c(-5e5,5e5))
plot(data$Distance.to.TSS, data$tagsPerInterval, type='h', xlab='Distance to transcriptional start site',
     ylab='Tags per base pair', main='Density of CREB binding', xlim=c(-5e5,5e5))
par(mfrow=c(2,1))
plot(stats::density(data$Distance.to.TSS,bw='nrd0'), xlim=c(-1e4, 1e4),main="CREB binding at TSS",xlab="Distance from TSS (upstream - downstream)")
polygon(density(data$Distance.to.TSS,bw='nrd0'), col="lightblue", border="grey") 
plot(stats::density(data$intervalSize,bw='nrd0'), xlim=c(0, 1e4),main="Proximal promoter length", xlab="Length of proximal promoter")
polygon(density(data$intervalSize,bw='nrd0'), col="lightgreen", border="grey") 
par(mfrow=c(1,1))

plot(data$CREBmotifNo,data$PeakScore,pch=20, col=rainbow(20),main='Tag count vs CREB motifs', xlab='Number of CREB motifs',ylab='Tag counts')

x <- rnorm(data$CREBmotifNo)
y <- rnorm(data$PeakScore)
bin<-hexbin(x, y, xbins=100)
plot(bin, main="Hexagonal Binning") 

genome <- BSgenome.Hsapiens.UCSC.hg19
len = as.vector(seqlengths(genome)[1:24])
bigRange = GRanges(seqnames=Rle(data$Chr), ranges=IRanges(start=data$Start, end=data$End,names=data$Gene.Name),
                   strand=data$Strand, peakScore=data$PeakScore, TSSdist=data$Distance.to.TSS, motifs=data$CREBmotifNo,
                   motifDensity=data$motifPerInterval, tagDensity=data$tagsPerInterval)


#write.table(dataMotifDensity, './130322_CREBchipAnnotationFiltered.txt', sep='\t', row.names=F)