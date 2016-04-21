#Forest plots of hazard ratio
library(rmeta)
data = read.delim('~/Documents/stemCellSig/130117_signature/130130_stemCellSigTCGA_Cox.txt', skip=3, row.names=1)
data = data[1:15,]
labels=read.delim('~/Documents/stemCellSig/130117_signature/130219_stemCellSigTCGA_Labels.txt')

#must input the yaxis labels as a matrix
labels=as.matrix(labels)

forestplot(labels, data$exp.coef., data$lower.095, data$upper.95,zero=1,xlab='Hazard ratio',align=c('l','r','c'),boxsize=0.3)