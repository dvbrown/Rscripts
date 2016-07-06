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

root <- "/Users/u0107775/"
writedir <- "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/Coverage_GC/"
window = "500000"
binsize <- "500K"
name = as.character(args[5])
readlength = 63
input <- "Single-cell"
type <- paste(readlength,"bases_mappable",sep="")

library("limma")

setwd(writedir)
#samplelist <- read.table("sample.txt", header=F, sep="\t")
samplelist = list.files("~/Data/GandT_Seq/160606_ValidationGandT/Coverage_Raw/")

setwd(root)
### read the GC content file (contains 1:10002-20000	0.5923592359235924)
gccontent <- read.table(paste("~/Bioinformatics/resources/Mappablebins/BINS_11Sep2015/","Combined_Human_NCBI37.2_",binsize,"_",type,"_bins_GCperc_INPUT.txt",sep=""),
                        header=F,sep="\t")
colnames(gccontent)<-c("interval","GC","abs")

setwd(root)
gc <- strsplit(as.character(gccontent$interval),split=":")
gc <- do.call("rbind",gc)
startstop <- strsplit(gc[,2],split="-")
startstop <- do.call("rbind", startstop)
rownames(gccontent) <- paste(as.character(gc[,1]),as.character(startstop[,2]),sep="_")
gc <- "NA"
startstop <-"NA"

#############################################################################################################
######### run all samples or the specific one
######### For a specific sample set item = FILENAME
#############################################################################################################

names(samplelist) = samplelist

for (item in samplelist) {
  sampletostudy = as.character(samplelist[item])
  
  setwd(writedir)
  # sample <- read.table(paste(sampletostudy,"_REF_",binsize,"_mappable_",readlength,"bases.sorted.count",sep=""), header=T, sep="\t")
  sample <- read.table(paste("~/Data/GandT_Seq/160606_ValidationGandT/Coverage_Raw/", sampletostudy, sep=""), header=F, sep="\t")
  colnames(sample) = c('chromosome', 'start', 'end', 'x', 'bin', 'test')
  
  
  ### give new rownames to sample
  rownames(sample) <- paste(as.character(sample$chromosome),as.character(sample$end),sep="_")
  ### intersect them
  same <- intersect(rownames(sample),rownames(gccontent))
  
  ###now ADD gccontent column to "sc"
  sc <- cbind(sample[same,],gccontent[same,2],gccontent[same,3])
  ### assign the new added column name as "GC"
  colnames(sc) <- c(colnames(sample),"GC","abs")
  
  ### filter on GC-content: removes regions with GC-contents lower than 28%
  scgcfiltered <- sc[as.numeric(sc$GC)>0.28,]
  nrow(scgcfiltered)
  
  ### perform GC-correction on raw log2ratio
  a <- as.numeric(scgcfiltered$test) + 1
  ratio <- a / mean(a)
  
  ### store log2 value of logR in logRprior
  logRprior <- log2(ratio)
  ### store the filtered file's GC content in gcplot
  gcplot <- scgcfiltered$GC
  
  
  ### apply loessfit to logRprior and gcplot with different spans and store them
  fitgenome1 <- loessFit(logRprior,gcplot,span=1)
  fitgenome05 <- loessFit(logRprior,gcplot,span=0.5)
  fitgenome <- loessFit(logRprior,gcplot)
  
  ### for GC correction substract "fitgenome$fitted" for respective logR... now again substract median of "logRprior-fitgenome$fitted)
  logR <- (logRprior-fitgenome$fitted)-median((logRprior-fitgenome$fitted)[as.character(scgcfiltered$chromosome)!="Y" & as.character(scgcfiltered$chromosome)!="X"],na.rm=T)
  logR[is.na(logR)] <- logRprior[is.na(logR)]
  
  # Plot
  setwd(writedir)
  jpeg(paste(sampletostudy,".loessFITonLOGRprior-",binsize,".",type,"-M30-II.jpg",sep=""),pointsize=12,width=1600,height=900)
  print(paste(sampletostudy,".loessFITonLOGRprior-",binsize,".",type,"-M30-II.jpg",sep=""))
  
  par(mfrow=c(1,2))
  plot(gcplot,logRprior,pch=19,cex=0.1,main=paste(sampletostudy,"Genome-wide GCcontent-influence",sep="_"), xlab="GC-content")
  points(gcplot,fitgenome1$fitted, pch=19, cex=0.5, col='orange')
  points(gcplot,fitgenome05$fitted, pch=19, cex=0.5, col='pink')
  points(gcplot,fitgenome$fitted, pch=19, cex=0.5, col='red')
  abline(h=0,col="grey")
  plot(gcplot,logR,pch=19,cex=0.1,main=paste(sampletostudy,"Genome-wide GC-corrected",sep="_"), xlab="GC-content")
  print(paste(sampletostudy,"Genome-wide GC-corrected",binsize,sep="_"))
  points(gcplot,fitgenome1$fitted-fitgenome$fitted, pch=19, cex=0.5, col='orange')
  points(gcplot,fitgenome05$fitted-fitgenome$fitted, pch=19, cex=0.5, col='pink')
  points(gcplot,fitgenome$fitted-fitgenome$fitted, pch=19, cex=0.5, col='red')
  abline(h=0,col="grey")
  dev.off()
  
  chr <- scgcfiltered$chr
  start <- scgcfiltered$start
  end <- scgcfiltered$end
  gcplot <- scgcfiltered$GC
  abs <- scgcfiltered$abs
  
  # Write output
  setwd(writedir)
  towritetable <- cbind(as.character(chr),start,logR,end,abs)
  colnames(towritetable) <- c("chromosome","position","log2","end","abs")
  print(paste(sampletostudy,".M30.hits-vs-NoREF.M30.hits.window-",window,".",type,".minw-4.GCcorrected-M30-II.cnv",sep=""))
  write.table(towritetable, paste(sampletostudy,".M30.hits-vs-NoREF.M30.hits.window-",window,".",type,".minw-4.GCcorrected-M30-II.cnv",sep=""), quote=F, col.names=T, row.names=F, sep="\t")
  
}  