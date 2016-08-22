cnv <- vector("list",300)
mapd <- vector("list",300)


thresholds <- paste("-PCF-kmin_3-gamma_",gamma,"-fastPCF-window-",window,"vs-NoREF-M30-II",sep="")
for(input in sampleSTART:sampleEND){
setwd(writedir)
readlength = as.numeric(as.character(samplelist[input,2]))
type <- paste(readlength,"bases_mappable",sep="")
print(paste("SEGMENTSlogR.GCcorrected-M30-",as.character(samplelist[input,1]),thresholds,type,".txt",sep=""))

###### Upload the LogR file
cnv[[input]] <- read.table(paste("SEGMENTSlogR.GCcorrected-M30-",as.character(samplelist[input,1]),thresholds,type,".txt",sep=""), sep="\t", header=T)

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


setwd(QualityCheck)
mappedplot <- do.call("rbind",mapped[sampleSTART:sampleEND])
mappedplot <- cbind(as.data.frame(samplelist[sampleSTART:sampleEND,1]),mappedplot)
colnames(mappedplot) <- c("name","sequenced","mapped","MappedPerc","Depth","Breadth", "MAPD")
write.table(mappedplot,paste("Quality_Check_Mapping_Stats_from_",sampleSTART,"to",sampleEND,"_",name,"_PairedEnd-vs-NoREF.gamma",gamma,".",binsize,".",type,".txt",sep=""),col.names=TRUE, row.names=FALSE, sep="\t",quote=F,dec=",")
}

