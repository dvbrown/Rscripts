##############################################################
#######  Steps to change
#######  1) root and writing dir 
#######  2) Edit "sample.txt".. i have attached the sample.txt as an example
#######  3) Edit Sample to study "sampletostudy"... which is on the basis of your "filenames"
#######  4) check whether the GC bins files and Mappable bins files are present
#######  Done
##############################################################
#######  Removed the chromosome Y while calculating the segmentation
###############################################################

root <- "Mention your working Directory"
writedir <- "Mention your wrinitng Directory"
window = "Mention window size e.g. 500000"
binsize <- "Mention window size e.g. 500K"
name = as.character(args[5])
readlength = 69
input <- "Mention the Cell number"
type <- paste(readlength,"bases_mappable",sep="")

library("limma")

setwd(writedir)
samplelist <- read.table("sample.txt", header=F, sep="\t")

library("limma")
source("/uz/data/avalok/symbiosys/raw/Parveendabas/bwaHumanIndexfile/OWN_CNV_count/fastPCF.R")

sampletostudy = as.character(samplelist[input,1])

setwd(writedir)
for(gamma in c(25))
{

    setwd(writedir)    
    thresholds <- paste("-PCF-kmin_3-gamma_",gamma,"-fastPCF-window-",window,"vs-NoREF-M30-II",sep="")
    
    scdata <- read.delim(paste(sampletostudy,".M30.hits-vs-NoREF.M30.hits.window-",window,".",type,".minw-4.GCcorrected-M30-II.cnv",sep=""))
    
    chrom <- c(seq(1,22),"X","Y")
    SEGMENTATIONtable <- vector("list",length(chrom))


    for(chr in 1:23)
    {
      print(paste(type," chr ",chr," of binsize ",binsize," ",sampletostudy,sep=""))
      pos <- scdata$position[as.character(scdata$chromosome)==chrom[chr]]
      posend <- scdata$end[as.character(scdata$chromosome)==chrom[chr]]
      posabs <- scdata$abs[as.character(scdata$chromosome)==chrom[chr]]
      data <- scdata$log2[as.character(scdata$chromosome)==chrom[chr]]
      y <- data[order(pos)][is.finite(data[order(pos)])]
      posselect <- pos[order(pos)][is.finite(data[order(pos)])]
      posendselect <- posend[order(posend)][is.finite(data[order(posend)])]
      posabsselect <- posabs[order(posabs)][is.finite(data[order(posabs)])]
      
      sdev=getMad(y,k=25)
      res=selectFastPcf(y,3,gamma*sdev,T)
      segments=res$yhat
      
      SEGMENTATIONtable[[chr]] <- cbind(rep(chrom[chr],length(posselect)), posselect, y, segments,posendselect,posabsselect)
    }

    towriteSEGMENTS <- do.call("rbind", SEGMENTATIONtable)
    colnames(towriteSEGMENTS) <- c("Chr","Pos","logR","logRsegment","end","abs")
    
    setwd(writedir)
    write.table(towriteSEGMENTS, paste("SEGMENTSlogR.GCcorrected-M30-",as.character(samplelist[input,1]),thresholds,type,".txt",sep=""), quote=F,sep="\t",col.names=T,row.names=F)

	#######################################
	### plotting logR SEGMENTATION
	setwd(writedir)
	for(chr in 1:23){

	jpeg(paste(sampletostudy,"Complete_chr_logR.gamma",gamma,".",type,binsize,"_",chr,".jpeg",sep=""),width = 1200, height = 600, units = "px", pointsize = 20)		#,width=10,height=10)
	#par(mfrow=c(2,1))
	limx <- max(as.numeric(as.character(towriteSEGMENTS[as.character(towriteSEGMENTS[,1])==chrom[chr],2]))[is.finite(as.numeric(as.character(towriteSEGMENTS[as.character(towriteSEGMENTS[,1])==chrom[chr],2])))])


	plot(as.numeric(as.character(towriteSEGMENTS[as.character(towriteSEGMENTS[,1])==chrom[chr],2])),as.numeric(as.character(towriteSEGMENTS[as.character(towriteSEGMENTS[,1])==chrom[chr],3])),col="black",pch=19,ylim=c(-4,4),xlim=c(1,limx),cex=0.5, xlab=NA,ylab=paste("logR (chr ",chrom[chr],") ",sep=""),main=paste(as.character(samplelist[input,1]),".",type,binsize,".gamma",gamma,sep=""),cex.main=1)
	points(as.numeric(as.character(towriteSEGMENTS[as.character(towriteSEGMENTS[,1])==chrom[chr],2])),as.numeric(as.character(towriteSEGMENTS[as.character(towriteSEGMENTS[,1])==chrom[chr],4])),col="yellow",pch=19,cex=0.7)
	dev.off()
	}


   }
