
### Load the LogR/CN file
setwd(LogRdir)
thresholds <- paste("-PCF-kmin_3-gamma_",gamma,"-fastPCF-window-",window,"vs-NoREF-M30-II",readlength,"bases_mappable",sep="")
print(paste("SEGMENTSlogR.GCcorrected-M30-",as.character(samplelist[input,1]),thresholds,".txt",sep=""))
cnv[[input]] <- read.table(paste("SEGMENTSlogR.GCcorrected-M30-",as.character(samplelist[input,1]),thresholds,".txt",sep=""), sep="\t", header=T)

### Get the chromosome boundaries
        chr <- as.character(cnv[[input]]$Chr)
        chr.shift <- c("chrom",chr[-(length(chr))])
        vlines <- c(1,cnv[[input]]$abs[which(chr != chr.shift) -1], cnv[[input]]$abs[nrow(cnv[[input]])])
        chr.text <- c(1:22, "X", "Y")
        chr.text <- chr.text[c(1:14, 16, 18, 20, 22:24)]
        vlines.shift <- c(vlines[-1], 4*10^9)
        chr.at <- vlines + (vlines.shift - vlines) / 2
        chr.at <- chr.at[c(1:14, 16, 18, 20, 22:24)]

### Plot LogR/CN Genomewide
setwd(Combined)
plot(cnv[[input]][,6],as.numeric(as.character(cnv[[input]][,3])),col="grey",pch=19,ylim=c(-4,4),xlim=c(1,max(cnv[[input]]$abs)),cex=0.55, xlab=NA,ylab=paste("LogR",sep=""),main=paste(as.character(samplelist[input,1]),".gamma",gamma,".",type,binsize,sep=""),cex.main=1,cex.lab=1.3, cex.axis=1, cex.main=1.1, cex.sub=1.1,mgp = c(2, 0.5, 0))
        lines(x=cnv[[input]]$abs, y=cnv[[input]]$logR, col="#CCCCCC", cex=0.5)
        points(x=cnv[[input]]$abs, y=cnv[[input]]$logRsegment, col="red", cex=0.6)
        abline(v=vlines)
        abline(h=0,col="black",lwd=1,lty=2)
        mtext(chr.text, at = chr.at,cex=0.7)

