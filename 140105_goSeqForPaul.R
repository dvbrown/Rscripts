library(goseq)
setwd('~/Documents/CREB/paulEnrichment/')
intData = read.delim('140205_tcgaGeneMean_CREBsites')
intData = intData[!duplicated(intData$ENSEMBL),]

# up.genes <- as.integer(intData$rowMean > 2 & intData$CREBcounts > 0)
# names(up.genes) <- intData$ENSEMBL
# table(up.genes)
# 
# down.genes <- as.integer(intData$rowMean < 2 & intData$CREBcounts > 0)
# names(down.genes) <-  intData$ENSEMBL
# table(down.genes)

all.genes <- as.integer(abs(intData$rowMean) > 2 & intData$CREBcounts > 0)
names(all.genes) <-  intData$ENSEMBL
table(all.genes)

# The null model should be flat for microarray data
# null_model_up <- nullp(up.genes, genome = "hg19", id = "ensGene")
# null_model_down <- nullp(down.genes, genome = "hg19", id = "ensGene")

null_model_all <- nullp(all.genes, genome = "hg19", id = "ensGene")
# Calculate enrichment of KEGG enrichment. Discard full GO categories now

# kegg_enrichmentUp <- goseq(null_model_up, genome = "hg19", id = "ensGene", 
#                           test.cats = c("GO:BP","GO:MF"), method = "Wallenius")
# kegg_enrichmentUp <- goseq(null_model_up, genome = "hg19", id = "ensGene", 
#                            test.cats = "KEGG", method = "Wallenius")
# kegg_enrichmentDown <- goseq(null_model_down, genome = "hg19", id = "ensGene", 
#                              test.cats = "GO:MF", method = "Wallenius")
# kegg_enrichmentDown <- goseq(null_model_down, genome = "hg19", id = "ensGene", 
#                              test.cats = "KEGG", method = "Wallenius")

kegg_enrichmentAll <- goseq(null_model_all, genome = "hg19", id = "ensGene", 
                             test.cats = "GO:BP", method = "Wallenius")

# kegg_enrichmentAllKEGG <- goseq(null_model_all, genome = "hg19", id = "ensGene", 
# There were no enriched KEGG categories anyway for the full dataset        test.cats = "KEGG", method = "Wallenius")
# Control for FDR
# kegg_fdrs_up <- p.adjust(kegg_enrichmentUp[, 2], method = "BH")
# kegg_enrichment_up <- cbind(kegg_enrichmentUp, FDR = kegg_fdrs_up)
# 
# kegg_fdrs_down <- p.adjust(kegg_enrichmentDown[, 2], method = "BH")
# kegg_enrichment_down <- cbind(kegg_enrichmentDown, FDR = kegg_fdrs_down)

kegg_fdrs_all <- p.adjust(kegg_enrichmentAll[, 2], method = "BH")
kegg_enrichment_all <- cbind(kegg_enrichmentAll, FDR = kegg_fdrs_all)

# kegg_fdrs_allK <- p.adjust(kegg_enrichmentAllK[, 2], method = "BH")
# kegg_enrichment_allK <- cbind(kegg_enrichmentAllK, FDR = kegg_fdrs_allK)
# No KEGG pathways aere significant if I take all the genes regardless of up or downregulation

# Now annotate IDs
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("go_id", "name_1006"), #filters = "refseq_dna",
                 values = enriched_pathwaysAll$pathway, mart = mart)

# Annotate the KEGG IDs
library(KEGGREST)
keggid2name <- keggList("pathway", "hsa")
head(keggid2name)

# names(keggid2name) <- sapply(names(keggid2name), substring, 9)
# kegg_enrichment <- cbind(kegg_enrichment_up, pathway = keggid2name[kegg_enrichment_up$category])
# enriched_pathwaysUp <- kegg_enrichment[kegg_enrichment_up$FDR < 0.1, ]
# dim(enriched_pathwaysUp)[1]
# head(enriched_pathwaysUp)
# #write.table(enriched_pathwaysUp, './goSeq/140105_enrichedKEGGpathwaysShortTerm.txt', sep='\t')
# enriched_pathwaysUp1 = merge.data.frame(enriched_pathwaysUp, results, by.x='category', by.y='go_id')
# write.table(enriched_pathwaysUp1, './goSeq/140105_enrichedGOtermsUpGenes.txt', sep='\t',row.names=F)
# 
# keggid2name <- keggList("pathway", "hsa")
# head(keggid2name)
# names(keggid2name) <- sapply(names(keggid2name), substring, 9)
# kegg_enrichment <- cbind(kegg_enrichment_down, pathway = keggid2name[kegg_enrichment_down$category])
# enriched_pathwaysDown <- kegg_enrichment[kegg_enrichment_down$FDR < 0.1, ]
# dim(enriched_pathwaysDown)[1]
# head(enriched_pathwaysDown)
# #write.table(enriched_pathwaysDown, './goSeq/140105_enrichedKEGGpathwaysLongTerm.txt', sep='\t')
# enriched_pathwaysDown1 = merge.data.frame(enriched_pathwaysDown, results, by.x='category', by.y='go_id')
# write.table(enriched_pathwaysDown1, './goSeq/140105_enrichedGOtermsDownGenes.txt', sep='\t', row.names=F)

keggid2name <- keggList("pathway", "hsa")
head(keggid2name)
names(keggid2name) <- sapply(names(keggid2name), substring, 9)
kegg_enrichment <- cbind(kegg_enrichment_all, pathway = keggid2name[kegg_enrichment_all$category])
enriched_pathwaysAll <- kegg_enrichment #[kegg_enrichment_all$FDR < 0.1, ] KEEP ALL THE PATHWAYS
dim(enriched_pathwaysAll)[1]
head(enriched_pathwaysAll)
#write.table(enriched_pathwaysDown, './goSeq/140105_enrichedKEGGpathwaysLongTerm.txt', sep='\t')
enriched_pathwaysAll1 = merge.data.frame(enriched_pathwaysAll, results, by.x='category', by.y='go_id')
write.table(enriched_pathwaysAll1, './goSeq/140105_enrichedGOtermsAllGenesGO_BP.txt', sep='\t', row.names=F)

# Now divide the numDE in cat by numInCat
enriched_pathwaysAll1$percentDEofTotal = enriched_pathwaysAll1$numDEInCat / enriched_pathwaysAll1$numInCat *100
plot(enriched_pathwaysAll1$percentDEofTotal, -log10(enriched_pathwaysAll1$FDR), xlim=c(0,100))

##################################### Make a ggplot object ##################################################
source('~/Documents/Rscripts/120704-sortDataFrame.R')
enriched_pathwaysAll1 = sort.dataframe(enriched_pathwaysAll1, 6, highFirst=FALSE)
library(ggplot2))
enriched_pathwaysAll1$threshold = as.factor(abs(enriched_pathwaysAll1$percentDEofTotal > 15 & enriched_pathwaysAll1$FDR < 0.05))

g = ggplot(data=enriched_pathwaysAll1, aes(x=percentDEofTotal, y=-log10(FDR), colour=threshold)) +
  geom_point(alpha=0.80, size=2) +
  opts(legend.position = "none", title=("CREB dependent pathways disrupted in GBM")
  ) +
  #xlim(c(-5, 5)) + ylim(c(0, 50)) +
  xlab("Proportion of CREB regulated disrupted GBM genes in pathway") + ylab("Significance")
g
#subset gene names for only significant genes
# dd_text = enriched_pathwaysAll1[(abs(enriched_pathwaysAll1$percentDEofTotal) > 20) & (enriched_pathwaysAll1$FDR < 0.01),]
dd_text = enriched_pathwaysAll1[c(19, 36, 47),]

#add text to volcano
g + geom_text(data = dd_text, aes(x=percentDEofTotal, y=-log10(FDR),
                                  label=name_1006, size=0.2), colour="black")