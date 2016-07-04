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

root_logR <- "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/DNA/LogR/"
root_CN <- "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/DNA/Copy_number/"
writedir <- "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/DNA/Copy_number/WholeGenome/"
window = "500000"
binsize <- "500K"
#name = as.character(args[5])
readlength = 63
#input <- 1
type <- paste(readlength,"bases_mappable",sep="")
plCN = 2
plchr =12
plstart = 45000000
plend = 130000000

samplelist <- list.files("/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/DNA/LogR/", pattern = '*.txt')
gamma<-10

## rest
for(input in 1:length(samplelist)){
print(input)
 sampletostudy = as.character(input)
    print(sampletostudy)
 thresholds <- paste("-PCF-kmin_3-gamma_",gamma,"-fastPCF-window-",window,"vs-NoREF-M30-II",sep="")
    chrom <- c(seq(1,22),"X","Y")
   
setwd(root_logR)
  #Define segmentation threshold data
LogR <-read.table(sampletostudy,sep="\t",header=T)
 
##Plots 
#pdf("genomewideProfiles.pdf",onefile=T)

setwd(writedir)
titl = paste(as.character(substr(sampletostudy,1,69)),".refLOCUS",plchr,".",type,binsize,".gamma",gamma,sep="")
jpeg(paste(titl,".jpeg",sep=""),
     width = 1200, height = 600, units = "px", pointsize = 20)		#,width=10,height=10)

### plotProfile
#par(mfcol=c(2,1))
data<-LogR$logR
min<- min(data)
max<-6
chrom<-as.vector(LogR$Chr)
chrom<-as.vector(LogR$Chr)
chrom[chrom=="X"]<-23
chrom<-as.numeric(chrom)

plot(c(1:length(data)),data, pch=19,cex=0.3,col="black",type="p", ylim=c(min-0.5,max+0.5), xlab = "Chromosomes", ylab = "logR ", main=titl ,xaxt = "n",cex.lab=1.2)
points(c(1:length(data)),LogR$logRsegment, pch=19,cex=0.3,type="p", ylim=c(min-0.5,max+0.5),xaxt = "n",cex.lab=1.2,col="yellow")
#tempSeg<-segmentData[segmentData$Sample==i,]
#segments(c(1,cumsum(tempSeg$nclone)[-nrow(tempSeg)]), tempSeg$values,cumsum(tempSeg$nclone), tempSeg$values,lwd=5,col="black")
axis(side = 1, at = as.integer(c(1:length(data))[1]), label = "", tick = FALSE)
#for(j in 1:length(c)){
#tempSeg<-segmentData[segmentData$Sample==m[j],]
#segments(c(1,cumsum(tempSeg$nclone)[-nrow(tempSeg)]), tempSeg$values,cumsum(tempSeg$nclone), tempSeg$values,lwd=5-(j*1.2), col = c[j])
#}
ax <- (cumsum(table(chrom)) + c(0, cumsum(table(chrom))[-length(cumsum(table(chrom)))]))/2 
axis(side = 1, at = ax, labels = c(1:23), cex = 0.2,lwd = 0.5, las = 1, cex.axis = 1, cex.lab = 1)
abline(h = 0)
for (iii in 1:length(cumsum(table(chrom)))) {segments(cumsum(table(chrom))[[iii]], -100, cumsum(table(chrom))[[iii]],100, lty = 2) }
graphics.off()

}


#Copy number
setwd(root_CN)
samplelist <- list.files("/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/DNA/Copy_number/", pattern = '*.txt')

for(input in samplelist){
  print(input)
  sampletostudy = as.character(input)
  print(sampletostudy)
  thresholds <- paste("-PCF-kmin_3-gamma_",gamma,"-fastPCF-window-",window,"vs-NoREF-M30-II",sep="")
  chrom <- c(seq(1,22),"X","Y")

  CN<-read.table(sampletostudy,sep="\t",header=T)
  
  jpeg(paste("../genomewideProfilesCN500K/Copynumber_", gsub(".sorted_dedup._mappable_500K_89bases.sorted.count","",sampletostudy),".jpeg",sep=""),width = 1200, height = 600, units = "px", pointsize = 20)	
  
  data<-CN$rawCN
  min<- min(data)
  max<-max(10)
  #chrom<-as.vector(data$Chr)
  title<-paste("Copy number", gsub(".sorted_dedup._mappable_500K_89bases.sorted.count","",sampletostudy))
  plot(c(1:length(data)),data, pch=19,cex=0.3,col="black",type="p", ylim=c(min-0.5,max+0.5), xlab = "Chromosomes", ylab = "Copy number ", main=title ,xaxt = "n",cex.lab=1.2)
  points(c(1:length(data)),CN$segmentedCN, pch=19,cex=0.3,type="p", ylim=c(min-0.5,max+0.5),xaxt = "n",cex.lab=1.2,col="yellow")
  #tempSeg<-segmentData[segmentData$Sample==i,]
  #segments(c(1,cumsum(tempSeg$nclone)[-nrow(tempSeg)]), tempSeg$values,cumsum(tempSeg$nclone), tempSeg$values,lwd=5,col="black")
  axis(side = 1, at = as.integer(c(1:length(data))[1]), label = "", tick = FALSE)
  #for(j in 1:length(c)){
  #tempSeg<-segmentData[segmentData$Sample==m[j],]
  #segments(c(1,cumsum(tempSeg$nclone)[-nrow(tempSeg)]), tempSeg$values,cumsum(tempSeg$nclone), tempSeg$values,lwd=5-(j*1.2), col = c[j])
  #}
  ax <- (cumsum(table(chrom)) + c(0, cumsum(table(chrom))[-length(cumsum(table(chrom)))]))/2 
  axis(side = 1, at = ax, labels = c(1:23), cex = 0.2,lwd = 0.5, las = 1, cex.axis = 1, cex.lab = 1)
  abline(h =2)
  for (iii in 1:length(cumsum(table(chrom)))) {segments(cumsum(table(chrom))[[iii]], -100, cumsum(table(chrom))[[iii]],100, lty = 2) }
  graphics.off()

}

