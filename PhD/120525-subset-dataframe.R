#remove non-DE probes from the column of interest
results.strip0 = subset.data.frame(results.data.frame, gliomaStem.vs.normalStem != 0)

#remove differentially expressed genes from the other columns. Continue for combinations
results.strip1 = subset.data.frame(results.strip0, gliomaStem.vs.normalBrain == 0)

#read in gene lst from probe annotation script
genelist.DataFrame = as.data.frame(genelist)

#I had to fix column names of results.strip in excel to ensure column header == 'GeneID'
merged.results = merge.data.frame(results.strip1, genelist.DataFrame)
head(merged.results)
write.table(merged.probe, './120524-reanalysis-output/120525-geneID-diffgenes-GSCvsNSConly.txt', sep = '\t', quote=F)