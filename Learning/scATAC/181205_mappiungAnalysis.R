setwd("~/Data/")

a2 = read.delim("mappingStats_a2.tsv",header=F)
a2$sample = "A2"
a2Unmapped = 185981
a2 = a2[c(1:286),]
a3 = read.delim("mappingStats_a3.tsv", header = F)
a3$sample = "A3"
a3Unmapped = 202125
a3 = a3[c(1:286),]

dat = rbind(a2,a3)
colnames(dat) = c("chr", "size", "reads", "unmap", "sample")

colnames(a2) = c("chr", "size", "reads", "unmap", "sample")
colnames(a3) = c("chr", "size", "reads", "unmap", "sample")

# Mitochondrial %
a2Mito = a2[25,3] / sum(a2$reads)
a3Mito = a3[25,3] / sum(a3$reads)

# ERCC %
a2ERCC = sum(a2[c(195:286),3]) / sum(a2$reads)
a3ERCC = sum(a3[c(195:286),3]) / sum(a3$reads)

# Unmapped %
a2UmapPercent = a2Unmapped / (sum(a2[,"reads"]) + a2Unmapped)
a3UmapPercent = a3Unmapped / (sum(a3[,"reads"]) + a3Unmapped)

# Make into a nice dataframe
output = as.data.frame(cbind(rbind(a2UmapPercent, a3UmapPercent),
                             rbind(a2Mito, a3Mito),
                             rbind(a2ERCC, a3ERCC)))
output = output *100
colnames(output) = c("Unmapped", "mtDNA%", "ERCC%")
row.names(output) = c("A2","A3")
write.table(output, "./181206_libraryStats.tsv", sep = "\t", quote = F)