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

root_CN <- "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/Copy_number/"
writedir <- "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/Copy_number/WholeGenomeCNV/"
window = "500000"
binsize <- "500K"
#name = as.character(args[5])
readlength = 63
type <- paste(readlength,"bases_mappable",sep="")
plCN = 3
plchr =12
plstart = 45000000
plend = 130000000

gamma<-10

####    Copy number    ####
setwd(root_CN)
samplelist <- list.files(pattern = '*.txt')

input = "SEGMENTSlogR.GCcorrected-M30-GC032370_AAGAGGCA-AGGCTTAG.cov.txt.M30.hits-vs-NoREF.M30.hits.window-500000.63bases_mappable.minw-4.GCcorrected-M30-II.cnv-PCF-kmin_3-gamma_25-fastPCF-window-500000vs-NoREF-M30-II75bases_mappable.txt"

for(input in samplelist){
  
  sampletostudy = as.character(input)
  print(input)
  thresholds <- paste("-PCF-kmin_3-gamma_",gamma,"-fastPCF-window-",window,"vs-NoREF-M30-II",sep="")
  chrom <- c(seq(1,22),"X")
  setwd(root_CN)
  CN<-read.table(input,sep="\t",header=T)
  titl = paste(as.character(substr(sampletostudy,1,16)), as.character(substr(sampletostudy,85,101)),
               ".refLOCUS",plchr,".",type,binsize,".gamma",gamma,sep="")
  
  setwd(writedir)
  jpeg(paste("copyNumber_", titl, ".jpeg",sep=""),
       width = 1200, height = 600, units = "px", pointsize = 20)		#,width=10,height=10)
  
  data<-CN$rawCN
  min<- min(data)
  max<-max(10)
  
  plot(c(1:length(data)),data, pch=19,cex=0.3,col="black",type="p", ylim=c(min-0.5,max+0.5), xlab = "Chromosomes", ylab = "Copy number ", main=substr(titl, 1, 33) ,xaxt = "n",cex.lab=1.2)
  points(c(1:length(data)),CN$segmentedCN, pch=19,cex=0.3,type="p", ylim=c(min-0.5,max+0.5),xaxt = "n",cex.lab=1.2,col="yellow")
  
  axis(side = 1, at = as.integer(c(1:length(data))[1]), label = "", tick = FALSE)

  ax <- (cumsum(table(chrom)) + c(0, cumsum(table(chrom))[-length(cumsum(table(chrom)))]))/2 
  axis(side = 1, at = ax, labels = c(1:23), cex = 0.2,lwd = 0.5, las = 1, cex.axis = 1, cex.lab = 1)
  abline(h =2)
  for (iii in 1:length(cumsum(table(chrom)))) {segments(cumsum(table(chrom))[[iii]], -100, cumsum(table(chrom))[[iii]],100, lty = 2) }
  
  graphics.off()
  
}