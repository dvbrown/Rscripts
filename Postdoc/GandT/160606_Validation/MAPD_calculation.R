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

root_logR <- "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/LogR/"
writedir <- "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/Copy_number/WholeGenome/"
QualityCheck = "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/LogR/QualityCheck"
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

cnv <- vector("list",300)
mapd <- vector("list",300)


thresholds <- paste("-PCF-kmin_3-gamma_",gamma,"-fastPCF-window-",window,"vs-NoREF-M30-II",sep="")
####    Calculate MAPD    ####
setwd(root_logR)
samplelist <- list.files(pattern = '*.txt')

# Intialise a vector to store the results
result = vector(mode="numeric", length = length(samplelist))
names(result) = samplelist

for(input in samplelist){
  
#input = samplelist[1]  

type <- paste(readlength,"bases_mappable",sep="")
print(paste("SEGMENTSlogR.GCcorrected-M30-",as.character(samplelist[input]),thresholds,type,".txt",sep=""))

###### Upload the LogR file
setwd(root_logR)
cnv[[input]] <- read.table(input, sep="\t", header=T)

dist <- vector("list", nrow(cnv[[input]])-1)
for(i in 1:length(dist))
{
#print(i)

		#### Calculate the Difference of consecutive logR values
        if(as.character(cnv[[input]][i+1,1])==as.character(cnv[[input]][i,1])){
        dist[[i]]<- abs(cnv[[input]]$logR[i+1]-cnv[[input]]$logR[i])
        }
}

### Take the median of MAPD and store it
mapd[[input]] <- median(do.call("rbind",dist), na.rm=T)
rm(dist)

# Write the mapd into the result vector
result[input] = mapd[[input]]

# setwd(QualityCheck)
# mappedplot <- do.call("rbind",mapped[sampleSTART:sampleEND])
# mappedplot <- cbind(as.data.frame(samplelist[sampleSTART:sampleEND,1]),mappedplot)
# colnames(mappedplot) <- c("name","sequenced","mapped","MappedPerc","Depth","Breadth", "MAPD")
# write.table(mappedplot,paste("Quality_Check_Mapping_Stats_from_",sampleSTART,"to",sampleEND,"_",name,
#                              "_PairedEnd-vs-NoREF.gamma",gamma,".",binsize,".",type,".txt",sep=""),
#             col.names=TRUE, row.names=FALSE, sep="\t",quote=F,dec=",")

}

# Write the result into a datframe and save it
setwd(QualityCheck)
result = result[c(1:27)]
resultDF = t(as.data.frame(result))
colnames(resultDF) = "MAPD"
write.table(resultDF, file="mapdScore.txt", row.names = T, sep="\t")