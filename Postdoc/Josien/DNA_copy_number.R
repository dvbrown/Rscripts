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

plCN = "Mention the ploidy"
plchr = "Mention the chromosome"
plstart = "Mention the start"
plend = "Mention the end"

library("limma")
source("/uz/data/avalok/symbiosys/raw/Parveendabas/bwaHumanIndexfile/OWN_CNV_count/fastPCF.R")

setwd(writedir)
samplelist <- read.table("sample.txt", header=F, sep="\t")


#############################################################################################################
######### run all samples or the specific one
#############################################################################################################
for(gamma in c(25))
{

    setwd(writedir)    
    thresholds <- paste("-PCF-kmin_3-gamma_",gamma,"-fastPCF-window-",window,"vs-NoREF-M30-II",sep="")
    chrom <- c(seq(1,22),"X","Y")
   
    sampletostudy = as.character(samplelist[input,1])
    print(sampletostudy)

  #Define segmentation threshold data
  LogR <- read.table(paste("SEGMENTSlogR.GCcorrected-M30-",as.character(samplelist[input,1]),thresholds,type,".txt",sep=""), sep="\t", header=T)
  
  # define ploidy based on a certain segment:


  LogRlev=mean(LogR[as.character(as.vector(LogR[,1]))==as.character(plchr)&LogR[,2]>plstart&LogR[,2]<plend,4])
  ploidy = plCN/2^LogRlev

  
  CN = matrix(ncol = 6, nrow = dim(LogR)[1])
  colnames(CN) = c("Chr","Pos","rawCN","segmentedCN","end","abs")
  rownames(CN) = rownames(LogR)
  CN[,1] = as.vector(LogR[,1])
  CN[,2] = LogR[,2]
  CN[,3] = 2^(LogR[,3]+log(ploidy,2))
  CN[,4] = pmax(round(2^(LogR[,4]+log(ploidy,2))),0)
# CN[,4] = 2^(LogR[,4]+log(ploidy,2))
  CN[,5] = LogR[,5]
  CN[,6] = LogR[,6]
  
  write.table(CN,paste("SEGMENTSlogR.GCcorrected-M30-",as.character(samplelist[input,1]),thresholds,type,".copynumber.refLOCUS.",plchr,".txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
  print(paste("SEGMENTSlogR.GCcorrected-M30-",as.character(samplelist[input,1]),thresholds,type,".copynumber.refLOCUS.",plchr,".txt",sep="")) 

	#######################################
	### plotting CN SEGMENTATION
	setwd(writedir)
	for(chr in 1:23){

	jpeg(paste(sampletostudy,"Complete_chr_CN.gamma",gamma,".refLOCUS.",plchr,"_",type,binsize,"_",chr,".jpeg",sep=""),width = 1200, height = 600, units = "px", pointsize = 20)		#,width=10,height=10)
	#par(mfrow=c(2,1))
	limx <- max(as.numeric(as.character(CN[as.character(CN[,1])==chrom[chr],2]))[is.finite(as.numeric(as.character(CN[as.character(CN[,1])==chrom[chr],2])))])

	plot(as.numeric(as.character(CN[as.character(CN[,1])==chrom[chr],2])),as.numeric(as.character(CN[as.character(CN[,1])==chrom[chr],3])),col="black",pch=19,ylim=c(0,10),xlim=c(1,limx),cex=0.5, xlab=NA,ylab=paste("DNA copy number (chr ",chrom[chr],") ",sep=""),main=paste(as.character(samplelist[input,1]),".refLOCUS.",plchr,".",type,binsize,".gamma",gamma,sep=""),cex.main=1)
	points(as.numeric(as.character(CN[as.character(CN[,1])==chrom[chr],2])),as.numeric(as.character(CN[as.character(CN[,1])==chrom[chr],4])),col="yellow",pch=19,cex=0.7)
	dev.off()
	}

 

   }

