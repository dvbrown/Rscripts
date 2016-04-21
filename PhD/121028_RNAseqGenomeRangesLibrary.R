setwd('~/Documents/RNAseqAnalysis/')
library(GenomicFeatures)
library(Rsamtools)

#load feature coordinates from Ensembl
txdb=makeTranscriptDbFromUCSC(genome='hg19',tablename='ensGene')
#extract the gene specific coordinates from ENsembl
tx_by_gene=transcriptsBy(txdb,'gene')
#load bam file
CD133p_B_Read=readBamGappedAlignments("s_4_ACAGTG_CD133p_B.bam.bamToFastq_R1.mergeBam.bam.sortSamtools.bam")
CD133n_B_Read=readBamGappedAlignments("s_4_CGATGT_CD133n_B.bam.bamToFastq_R1.mergeBam.bam.sortSamtools.bam")
CD133p_A_Read=readBamGappedAlignments("s_4_TGACCA_CD133p_A.bam.bamToFastq_R1.mergeBam.bam.sortSamtools.bam")
CD133p_C_Read=readBamGappedAlignments("s_4_GCCAAT_CD133p_C.bam.bamToFastq_R1.mergeBam.bam.sortSamtools.bam")
CD133n_A_Read=readBamGappedAlignments("s_4_ATCACG_CD133n_A.bam.bamToFastq_R1.mergeBam.bam.sortSamtools.bam")
CD133n_C_Read=readBamGappedAlignments("s_4_TTAGGC_CD133n_C.bam.bamToFastq_R1.mergeBam.bam.sortSamtools.bam")
CD133x_B_Read=readBamGappedAlignments("s_4_ACTTGA_CD133x_B.bam.bamToFastq_R1.mergeBam.bam.sortSamtools.bam")
CD133x_A_Read=readBamGappedAlignments("s_4_CAGATC_CD133x_A.bam.bamToFastq_R1.mergeBam.bam.sortSamtools.bam")
CD133x_C_Read=readBamGappedAlignments("s_4_GATCAG_CD133x_C.bam.bamToFastq_R1.mergeBam.bam.sortSamtools.bam")
CD133t_A_Read=readBamGappedAlignments("s_4_GGCTAC_CD133t_A.bam.bamToFastq_R1.mergeBam.bam.sortSamtools.bam")
CD133t_B_Read=readBamGappedAlignments("s_4_CTTGTA_CD133t_B.bam.bamToFastq_R1.mergeBam.bam.sortSamtools.bam")

#Check that the chromosome naming convention is the same for the reads and genome coordinates otherwise no matches
#The annotations have chromosomes called
names(seqlengths(tx_by_gene))
#The reads have chromosomes called
as.character(unique(rname(CD133nRead)))

countreads = function(GappedAlignments) { #Take a gapped alignments object (Bam file) and return a vector of counts.
  #change the bam file chromosome names to match the genome annotation, ie the tx_by_gene object
  new_read_chr_names=gsub("^","chr",rname(GappedAlignments))
  reads=GRanges(seqnames=new_read_chr_names,ranges=IRanges(start=start(GappedAlignments),end=end(GappedAlignments)), strand=rep("*",length(GappedAlignments)))
  #extract read counts
  counts=countOverlaps(tx_by_gene,reads)
  return (counts)
}

#Get the gene expression matrix for each sample
CD133p_B_count = countreads(CD133p_B_Read)
CD133n_B_count = countreads(CD133n_B_Read)
CD133p_A_count = countreads(CD133p_A_Read)
CD133p_C_count = countreads(CD133p_C_Read)
CD133n_A_count = countreads(CD133n_A_Read)
CD133n_C_count = countreads(CD133n_C_Read)
CD133x_B_count = countreads(CD133x_B_Read)
CD133x_A_count = countreads(CD133x_A_Read)
CD133x_C_count = countreads(CD133x_C_Read)
CD133t_A_count = countreads(CD133t_A_Read)
CD133t_B_count = countreads(CD133t_B_Read)

#encapsulate the separate read counts into a single dataframe
toc=data.frame(CD133n_rep1=CD133n_A_count, CD133n_rep2=CD133n_B_count, CD133n_rep3=CD133n_C_count,
               CD133p_rep1=CD133p_A_count, CD133p_rep2=CD133p_B_count, CD133p_rep3=CD133p_C_count,
               CD133x_rep1=CD133x_A_count, CD133x_rep2=CD133x_B_count, CD133x_rep3=CD133x_C_count,
               CD133t_rep1=CD133t_A_count, CD133t_rep2=CD133t_B_count,
               stringsAsFactors=FALSE)
rownames(toc)=names(tx_by_gene)
write.table(toc, './121028_readCountsGenomeRanges.txt', sep='\t')