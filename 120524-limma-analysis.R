setwd("~/Documents/public-datasets/microarray/GSC-NSC/Pollard2009")
data = read.table(file='110521-rmanormalize-data.txt', header=T)
colnames(data)

#read targets file
setwd('Rfiles/')
targets = readTargets('targets.txt', sep='\t')
f = paste(targets$source)
f = factor(f)
design = model.matrix(~0+f)
colnames(design) = levels(f)

fit = lmFit(data, design)
names(fit)

cont.matrix = makeContrasts(gliomaStem.vs.normalStem="glioma.stem-normal.stem", gliomaStem.vs.normalBrain=
  'glioma.stem-normal.brain', 
  normalStem.vs.normalBrain='normal.stem-normal.brain', levels=design)

fit2  = contrasts.fit(fit, cont.matrix)
fit2  = eBayes(fit2)

#write the output to a table of differentially expressed genes. Change this value to suit
gliomaStem.vs.normalStem.anno = topTable(fit2, coef=1, number=300, genelist=genelist,  
                                           adjust='BH', sort.by='logFC', lfc=1)

results <- decideTests(fit2)
summary(results)

#write the data out to a tab delimited text file
write.table(gliomaStem.vs.normalStem.anno, './output/110524-gliomastem-normalstem.txt', sep='\t', quote=F)

#write out some plots
par(mfrow=c(2,2))
vennDiagram(results, main='Natsume unpub upregulated genes', include='up', cex=0.6)
vennDiagram(results, main='Natsume unpub downregulated genes', include='down', cex=0.6)
vennDiagram(results, main='Natsume unpub total diff exoressed genes', include='both', cex=0.6)

volcanoplot(fit2, coef=1, highlight=10, names=annotaton$GeneSymbol, xlab='Log2 fold change', ylab='Negative log10 p-value', main='Glioma stem cells vs normal stem cells')