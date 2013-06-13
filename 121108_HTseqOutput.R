setwd('~/Documents/RNAdata/RNAseqAnalysis/121105_trimmomaticReads/mergedBam/121107_mergeSortTopHatAlignIndex/sortReadName/htSeqOut/')

CD133p_A_count = read.delim('s_4_TGACCA_CD133p_A_R1.fqpair.mergeBam.bam.sortSamtools.bam.sortReadName.txt', header=F, row.names=1)
CD133p_B_count = read.delim('s_4_ACAGTG_CD133p_B_R1.fqpair.mergeBam.bam.sortSamtools.bam.sortReadName.txt', header=F, row.names=1)
CD133p_C_count = read.delim('s_4_GCCAAT_CD133p_C_R1.fqpair.mergeBam.bam.sortSamtools.bam.sortReadName.txt', header=F, row.names=1)
CD133n_A_count = read.delim('s_4_ATCACG_CD133n_A_R1.fqpair.mergeBam.bam.sortSamtools.bam.sortReadName.txt', header=F, row.names=1)
CD133n_B_count = read.delim('s_4_CGATGT_CD133n_B_R1.fqpair.mergeBam.bam.sortSamtools.bam.sortReadName.txt', header=F, row.names=1)
CD133n_C_count = read.delim('s_4_TTAGGC_CD133n_C_R1.fqpair.mergeBam.bam.sortSamtools.bam.sortReadName.txt', header=F, row.names=1)
CD133x_B_count = read.delim('s_4_ACTTGA_CD133x_B_R1.fqpair.mergeBam.bam.sortSamtools.bam.sortReadName.txt', header=F, row.names=1)
CD133x_A_count = read.delim('s_4_CAGATC_CD133x_A_R1.fqpair.mergeBam.bam.sortSamtools.bam.sortReadName.txt', header=F, row.names=1)
CD133x_C_count = read.delim('s_4_GATCAG_CD133x_C_R1.fqpair.mergeBam.bam.sortSamtools.bam.sortReadName.txt', header=F, row.names=1)
CD133t_A_count = read.delim('s_4_GGCTAC_CD133t_A_R1.fqpair.mergeBam.bam.sortSamtools.bam.sortReadName.txt', header=F, row.names=1)
CD133t_B_count = read.delim('s_4_CTTGTA_CD133t_B_R1.fqpair.mergeBam.bam.sortSamtools.bam.sortReadName.txt', header=F, row.names=1)

toc=data.frame(CD133n_rep1=CD133n_A_count, CD133n_rep2=CD133n_B_count, CD133n_rep3=CD133n_C_count,
               CD133p_rep1=CD133p_A_count, CD133p_rep2=CD133p_B_count, CD133p_rep3=CD133p_C_count,
               CD133x_rep1=CD133x_A_count, CD133x_rep2=CD133x_B_count, CD133x_rep3=CD133x_C_count,
               CD133t_rep1=CD133t_A_count, CD133t_rep2=CD133t_B_count,
               stringsAsFactors=FALSE)

colnames(toc) = c('CD133n_rep1', 'CD133n_rep2', 'CD133n_rep3','CD133p_rep1', 'CD133p_rep2', 'CD133p_rep3', 'CD133x_rep1', 'CD133x_rep2','CD133x_rep3','CD133t_rep1', 'CD133t_rep2')