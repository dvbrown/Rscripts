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

root <- "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/DNA/LogR/"
writedir <- "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/Copy_number/"
window = "500000"
binsize <- "500K"
name = as.character(args[5])
readlength = 63
input <- "Single-cells"
type <- paste(readlength,"bases_mappable",sep="")

plCN = 2
plchr = '12'
plstart = '1'
plend = '60000'

library("limma")
source("~/Code/Rscripts/Postdoc/Josien/fastPCF.R")

setwd(writedir)
samplelist <- list.files(path="~/Data/GandT_Seq/160606_ValidationGandT/DNA/LogR/", pattern = ".txt")

#############################################################################################################
######### run all samples or the specific one
#############################################################################################################
for(gamma in c(25))
{

    setwd(root)   
    thresholds <- paste("-PCF-kmin_3-gamma_",gamma,"-fastPCF-window-",window,"vs-NoREF-M30-II",sep="")
    chrom <- c(seq(1,22),"X","Y")
   
    sampletostudy = as.character(samplelist[11])
    basename = substr(sampletostudy, 30, 56)
    print(sampletostudy)

  #Define segmentation threshold data
  LogR <- read.table(as.character(sampletostudy), sep="\t", header=T)
  
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
  
  setwd(writedir)
  write.table(CN,paste("SEGMENTSlogR.GCcorrected-M30-",as.character(basename),thresholds,type,".copynumber.refLOCUS.",plchr,".txt",sep=""),
              sep="\t",col.names=T,row.names=F,quote=F)
  print(paste("SEGMENTSlogR.GCcorrected-M30-",as.character(basename),thresholds,type,".copynumber.refLOCUS.",plchr,".txt",sep="")) 

	#######################################
	### plotting CN SEGMENTATION
	setwd(writedir)
	for(chr in 1:23){

	jpeg(paste(basename,"Complete_chr_CN.gamma",gamma,".refLOCUS.",plchr,"_",type,binsize,"_",chr,".jpeg",sep=""),
	     width = 1200, height = 600, units = "px", pointsize = 20)		#,width=10,height=10)
	#par(mfrow=c(2,1))
	limx <- max(as.numeric(as.character(CN[as.character(CN[,1])==chrom[chr],2]))[is.finite(as.numeric(as.character(CN[as.character(CN[,1])==chrom[chr],2])))])

	plot(as.numeric(as.character(CN[as.character(CN[,1])==chrom[chr],2])),as.numeric(as.character(CN[as.character(CN[,1])==chrom[chr],3])),
	     col="black",pch=19,ylim=c(0,10),xlim=c(1,limx),cex=0.5, xlab=NA,ylab=paste("DNA copy number (chr ",chrom[chr],") ",sep=""),
	     main=paste(as.character(samplelist[1]),".refLOCUS.",plchr,".",type,binsize,".gamma",gamma,sep=""),cex.main=1)
	
	points(as.numeric(as.character(CN[as.character(CN[,1])==chrom[chr],2])),
	       as.numeric(as.character(CN[as.character(CN[,1])==chrom[chr],4])),col="yellow",pch=19,cex=0.7)
	dev.off()
	}
   }




####    WHOLE GENOME PLOT ####
for (item in samplelist) {
  
  for(gamma in c(25))
  {
    
    setwd(root)   
    thresholds <- paste("-PCF-kmin_3-gamma_",gamma,"-fastPCF-window-",window,"vs-NoREF-M30-II",sep="")
    chrom <- c(seq(1,22),"X","Y")
    
    sampletostudy = item
    basename = substr(sampletostudy, 1, 69)
    print(sampletostudy)
    
    #Define segmentation threshold data
    LogR <- read.table(as.character(sampletostudy), sep="\t", header=T)
    
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
    
    setwd(writedir)
    write.table(CN,paste("SEGMENTSlogR.GCcorrected-M30-",as.character(basename),thresholds,type,".copynumber.refLOCUS.",plchr,".txt",sep=""),
                sep="\t",col.names=T,row.names=F,quote=F)
    print(paste("SEGMENTSlogR.GCcorrected-M30-",as.character(basename),thresholds,type,".copynumber.refLOCUS.",plchr,".txt",sep=""))
  }
    
  sampleName = item
  title = paste(as.character(substr(sampleName,1,69)),".refLOCUS",plchr,".",type,binsize,".gamma",gamma,sep="")
  
  jpeg(paste("WholeGenome/ploidy3/",title,".jpeg",sep=""),
       width = 1200, height = 600, units = "px", pointsize = 20)		#,width=10,height=10)
  #par(mfrow=c(2,1))
  limx <- max(as.numeric(as.character(CN[as.character(CN[,1])==chrom,2]))[is.finite(as.numeric(as.character(CN[as.character(CN[,1])==chrom,2])))])
  
  plot(as.numeric(as.character(CN[as.character(CN[,1])==chrom,2])),as.numeric(as.character(CN[as.character(CN[,1])==chrom,3])),
       col="black",pch=19,ylim=c(0,10),xlim=c(1,limx),cex=0.5, 
       xlab=CN[,1], # This needs to be fixed
       ylab="Copy Number",
       main=title,cex.main=1)
  
  points(as.numeric(as.character(CN[as.character(CN[,1])==chrom,2])),
         as.numeric(as.character(CN[as.character(CN[,1])==chrom,4])),col="yellow",pch=19,cex=0.7)
  dev.off()
}