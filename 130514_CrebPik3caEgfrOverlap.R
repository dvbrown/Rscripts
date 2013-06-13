library(VennDiagram)
library(venneuler)


currentDir = getwd()
setwd('~/Documents/CREB/targetGeneOverlap/DAVID/')
#obtained the list below from DAVID target gene list. Mainly text mining.
data = read.delim('130514_relatedGenesDAVID.txt', stringsAsFactors=FALSE, na.strings="")

l = list(creb=data$creb1, egfr=data$egfr, pik3ca=data$pik3ca)
creb = data$creb1
egfr = data$egfr
pi3k = data$pik3ca

venn.diagram(list(creb, egfr, pi3k), 
            filename='venn.tif', fill=c('red','blue', 'yellow'),main='This is a venn', euler.d=FALSE,
             na='remove')


#give up on vennDiagram. Only seems to return 1
options(java.parameters='-Xmx1192m')
#the intersect is the elements are in common = AND
#the union is elements in either list = OR
creb.egfr = intersect(data$creb1, data$egfr)
creb.pi3k = intersect(data$creb1, data$pik3ca)
egrf.pik3k = intersect(data$egfr, data$pik3ca)
overlap = intersect(creb.egfr, creb.pi3k)
A = length(data$creb1)
B = length(data$egfr)
C = length(data$pik3ca)

AB=length(creb.egfr)
AC=length(creb.pi3k)
BC=length(egrf.pik3k)
ABC=length(overlap)

vd <- venneuler(c(A=A, B=B, C=data$pik3ca, "A&B"=AB, 
                  "A&C"=AC, "B&C"=BC ,"A&B&C"=ABC))

vd <- venneuler(c(CREB1=97, EGFR=117, PI3K=80, "CREB1&EGFR"=70, "CREB1&PI3K"=52, "EGFR&PI3K"=59 ,"CREB1&EGFR&PI3K"=40))
plot(vd, main='Signaling pathway overlap\n DAVID downstream targets')

vd <- venneuler(c(CREB1=97, EGFR=117, PI3K=80, "CREB1&EGFR"=70, "CREB1&PI3K"=52, "EGFR&PI3K"=59 ,"CREB1&EGFR&PI3K"=40))
plot(vd, main='Signaling pathway overlap\n DAVID downstream targets')

write.table(overlap, '130515_overlap.txt', sep='\t', row.names=FALSE)