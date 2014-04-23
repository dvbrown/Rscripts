

############################################## Read in the resazurin assay readings ###############################################
setwd('~/Documents/Cell_biology/proliferation/Resazurin/140417_6clones/analysis/')
growthD3 = read.delim('140414_day3_linearRep.txt')
growthD7 = read.delim('140417_day7_linearRep.txt')

par(mfrow=c(2,1))
plot(growthD7$rep1, growthD7$rep2, ylab='replicate 2', xlab='replicate1', main='consistency day7')
plot(growthD3$rep1, growthD3$rep2, ylab='replicate 2', xlab='replicate1', main='consistency day3')
####################################################################################################################################



############################################## Read in the invasion assay readings #################################################
setwd('~/Documents/Cell_biology/microscopy/invasion/140414_invasion/')
invD3 = read.delim('140422_outputDay3.rep.txt')
invD7 = read.delim('140422_outputDay7.rep.txt')
####################################################################################################################################