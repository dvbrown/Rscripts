LogRdir = "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/DNA/LogR/"
Combined = "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/DNA/Copy_number/WholeGenome/"
copyNoDir = "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/DNA/Copy_number/"
samplelist <- list.files("/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/DNA/LogR/", pattern = '*.txt')

input = samplelist[1]

### Load the LogR/CN file
setwd(LogRdir)
thresholds <- paste("-PCF-kmin_3-gamma_",gamma,"-fastPCF-window-",window,"vs-NoREF-M30-II",readlength,"bases_mappable",sep="")
print(paste("processing",as.character(samplelist[1]),sep="__"))
cnv <- read.table(input, sep="\t", header=T)

### Get the chromosome boundaries
        chr <- as.character(cnv$Chr)
        chr.shift <- c("chrom",chr[-(length(chr))])
        vlines <- c(1,cnv$abs[which(chr != chr.shift) -1], cnv$abs[nrow(cnv)])
        chr.text <- c(1:22, "X", "Y")
        chr.text <- chr.text[c(1:14, 16, 18, 20, 22:24)]
        vlines.shift <- c(vlines[-1], 4*10^9)
        chr.at <- vlines + (vlines.shift - vlines) / 2
        chr.at <- chr.at[c(1:14, 16, 18, 20, 22:24)]

### Plot LogR/CN Genomewide ####
setwd(Combined)
sampleName = substr(input, 1, 69)

jpeg(paste(sampleName, "LogR",".jpeg",sep=""),width = 1200, height = 600, units = "px", pointsize = 20)	
plot(cnv[,6], as.numeric(as.character(cnv[,3])),
     col="grey",pch=19,ylim=c(-4,4),xlim=c(1,max(cnv$abs)),
     cex=0.55, xlab=NA,ylab=paste("LogR",sep=""),main=paste(as.character(sampleName),".gamma",gamma,".",type,binsize,sep=""),
     cex.main=1,cex.lab=1.3, cex.axis=1, cex.main=1.1, cex.sub=1.1,mgp = c(2, 0.5, 0))
        lines(x=cnv$abs, y=cnv$logR, col="#CCCCCC", cex=0.5)
        points(x=cnv$abs, y=cnv$logRsegment, col="red", cex=0.6)
        abline(v=vlines)
        abline(h=0,col="black",lwd=1,lty=2)
        mtext(chr.text, at = chr.at,cex=0.7)
dev.off()

#############################################
setwd(copyNoDir)
samplelist <- list.files("/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/DNA/Copy_number/", pattern = '*.txt')
input = samplelist[1]
cnv <- read.table(input, sep="\t", header=T)

### Get the chromosome boundaries
chr <- as.character(cnv$Chr)
chr.shift <- c("chrom",chr[-(length(chr))])
vlines <- c(1,cnv$abs[which(chr != chr.shift) -1], cnv$abs[nrow(cnv)])
chr.text <- c(1:22, "X", "Y")
chr.text <- chr.text[c(1:14, 16, 18, 20, 22:24)]
vlines.shift <- c(vlines[-1], 4*10^9)
chr.at <- vlines + (vlines.shift - vlines) / 2
chr.at <- chr.at[c(1:14, 16, 18, 20, 22:24)]

### Plot CN Genomewide ####
setwd(Combined)
sampleName = substr(input, 1, 69)

jpeg(paste(sampleName, "CopyNumber",".jpeg",sep=""),width = 1200, height = 600, units = "px", pointsize = 20)	
plot(cnv[,6], as.numeric(as.character(cnv[,4])),
     col="grey",pch=19,ylim=c(-4,4),xlim=c(1,max(cnv$abs)),
     cex=0.55, xlab=NA,ylab=paste("LogR",sep=""),main=paste(as.character(sampleName),".gamma",gamma,".",type,binsize,sep=""),
     cex.main=1,cex.lab=1.3, cex.axis=1, cex.main=1.1, cex.sub=1.1,mgp = c(2, 0.5, 0))
lines(x=cnv$abs, y=cnv$logR, col="#CCCCCC", cex=0.5)
points(x=cnv$abs, y=cnv$logRsegment, col="red", cex=0.6)
abline(v=vlines)
abline(h=0,col="black",lwd=1,lty=2)
mtext(chr.text, at = chr.at,cex=0.7)
dev.off()