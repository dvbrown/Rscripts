#Pool the 3 lists of CREB genes (gene expression, in silico CRE prediction and ChIP)
library(RankAggreg)
setwd('~/Documents/CREB/Zhang2005-CREBMontminyDB/')

expression = read.delim('GSE2060_RAW/human/HEK/120710-results/120711-forskolinVSaCREBforskolin.txt')

ChIP = read.delim('TableS4-ChIPChIP.txt')

CRE = read.delim('TableS1a-PromPredict.txt')

#trim data to interesting columns and have the same row number
CRE.ed = CRE[,c(1,2,3,6,7,10,13,15)]
ChIP.ed = ChIP[,c(2,3,5,8,9,10,11)]
expr.ed = expression[,c(1,2,3,4,6,7,8,9)]

#remove NAs and keep ony unique symbols
expr.unique = na.omit(expr.unique)
CRE.unique = na.omit(CRE.unique)
ChIP.unique = na.omit(ChIP.unique)

CRE.unique = unique(CRE.ed$Symbol)
ChIP.unique = unique(ChIP.ed$Symbol)
expr.unique = unique(expr.ed$GeneSymbol)

#extract the genesymbols, bind as rows
RankListInput = function (x,y,z) {
  list = data.frame(x, y, z)
  list1 = as.matrix(list)
  list.input = t(list1)
  list.input
                    }