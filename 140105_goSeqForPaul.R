library(goseq)
setwd('~/Documents/CREB/paulEnrichment/')
intData = read.delim('140205_tcgaGeneMean_CREBsites')
intData = intData[!duplicated(intData$ENSEMBL),]

up.genes <- as.integer(intData$rowMean > 2 & intData$CREBcounts > 0)
names(up.genes) <- intData$ENSEMBL
table(up.genes)

down.genes <- as.integer(intData$rowMean < 2 & intData$CREBcounts > 0)
names(down.genes) <-  intData$ENSEMBL
table(down.genes)

all.genes <- as.integer(abs(intData$rowMean) > 2 & intData$CREBcounts > 0)
names(all.genes) <-  intData$ENSEMBL
table(all.genes)

# The null model should be flat for microarray data
null_model_up <- nullp(up.genes, genome = "hg19", id = "ensGene")
null_model_down <- nullp(down.genes, genome = "hg19", id = "ensGene")

null_model_all <- nullp(all.genes, genome = "hg19", id = "ensGene")
# Calculate enrichment of KEGG enrichment. Discard full GO categories now

kegg_enrichmentUp <- goseq(null_model_up, genome = "hg19", id = "ensGene", 
                          test.cats = c("GO:BP","GO:MF"), method = "Wallenius")
#kegg_enrichmentUp <- goseq(null_model_up, genome = "hg19", id = "ensGene", 
#                           test.cats = "KEGG", method = "Wallenius")
kegg_enrichmentDown <- goseq(null_model_down, genome = "hg19", id = "ensGene", 
                             test.cats = c("GO:BP","GO:MF"), method = "Wallenius")
#kegg_enrichmentDown <- goseq(null_model_down, genome = "hg19", id = "ensGene", 
#                             test.cats = "KEGG", method = "Wallenius")

kegg_enrichmentAll <- goseq(null_model_all, genome = "hg19", id = "ensGene", 
                             test.cats = c("GO:BP","GO:MF"), method = "Wallenius")
kegg_enrichmentAllK <- goseq(null_model_all, genome = "hg19", id = "ensGene", 
                                                         test.cats = "KEGG", method = "Wallenius")
# Control for FDR
kegg_fdrs_up <- p.adjust(kegg_enrichmentUp[, 2], method = "BH")
kegg_enrichment_up <- cbind(kegg_enrichmentUp, FDR = kegg_fdrs_up)

kegg_fdrs_down <- p.adjust(kegg_enrichmentDown[, 2], method = "BH")
kegg_enrichment_down <- cbind(kegg_enrichmentDown, FDR = kegg_fdrs_down)

kegg_fdrs_all <- p.adjust(kegg_enrichmentAll[, 2], method = "BH")
kegg_enrichment_all <- cbind(kegg_enrichmentAll, FDR = kegg_fdrs_all)

kegg_fdrs_allK <- p.adjust(kegg_enrichmentAllK[, 2], method = "BH")
kegg_enrichment_allK <- cbind(kegg_enrichmentAllK, FDR = kegg_fdrs_allK)
# No KEGG pathways aere significant if I take all the genes regardless of up or downregulation

# Annotate the KEGG IDs
library(KEGGREST)
keggid2name <- keggList("pathway", "hsa")
head(keggid2name)

names(keggid2name) <- sapply(names(keggid2name), substring, 9)
kegg_enrichment <- cbind(kegg_enrichment_up, pathway = keggid2name[kegg_enrichment_up$category])
enriched_pathwaysUp <- kegg_enrichment[kegg_enrichment_up$FDR < 0.1, ]
dim(enriched_pathwaysUp)[1]
head(enriched_pathwaysUp)
#write.table(enriched_pathwaysUp, './goSeq/140105_enrichedKEGGpathwaysShortTerm.txt', sep='\t')
enriched_pathwaysUp1 = merge.data.frame(enriched_pathwaysUp, results, by.x='category', by.y='go_id')
write.table(enriched_pathwaysUp1, './goSeq/140105_enrichedGOtermsUpGenes.txt', sep='\t',row.names=F)

keggid2name <- keggList("pathway", "hsa")
head(keggid2name)
names(keggid2name) <- sapply(names(keggid2name), substring, 9)
kegg_enrichment <- cbind(kegg_enrichment_down, pathway = keggid2name[kegg_enrichment_down$category])
enriched_pathwaysDown <- kegg_enrichment[kegg_enrichment_down$FDR < 0.1, ]
dim(enriched_pathwaysDown)[1]
head(enriched_pathwaysDown)
#write.table(enriched_pathwaysDown, './goSeq/140105_enrichedKEGGpathwaysLongTerm.txt', sep='\t')
enriched_pathwaysDown1 = merge.data.frame(enriched_pathwaysDown, results, by.x='category', by.y='go_id')
write.table(enriched_pathwaysDown1, './goSeq/140105_enrichedGOtermsDownGenes.txt', sep='\t', row.names=F)

keggid2name <- keggList("pathway", "hsa")
head(keggid2name)
names(keggid2name) <- sapply(names(keggid2name), substring, 9)
kegg_enrichment <- cbind(kegg_enrichment_all, pathway = keggid2name[kegg_enrichment_all$category])
enriched_pathwaysAll <- kegg_enrichment[kegg_enrichment_all$FDR < 0.1, ]
dim(enriched_pathwaysAll)[1]
head(enriched_pathwaysAll)
#write.table(enriched_pathwaysDown, './goSeq/140105_enrichedKEGGpathwaysLongTerm.txt', sep='\t')
enriched_pathwaysAll1 = merge.data.frame(enriched_pathwaysAll, results, by.x='category', by.y='go_id')
write.table(enriched_pathwaysAll1, './goSeq/140105_enrichedGOtermsAllGenes.txt', sep='\t', row.names=F)

# Now annotate IDs
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("go_id", "name_1006"), #filters = "refseq_dna",
                 values = enriched_pathwaysUp$pathway, mart = mart)