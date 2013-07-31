library(MatrixEQTL)
setwd('~/Documents/eQTL/Matrix_eQTL_R/')
#This tutorial is hosted at http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/manual.html

#Matrix eQTL requires two input files, one with genotype information, and another with the gene expression. 
#It can also accept extra file with covariates and an error covariance matrix.

#The model of choice, there is more than 1 to choose from
useModel = modelLINEAR

#set up the SNP file and the expression file. The genotypes are coded as 0, 1 or 2 based on AA, AB, BB genotype.
SNP_file_name = 'Sample_Data/GE.txt'
expression_file_name = 'Sample_Data/SNP.txt'

#The covariate input file. This will require much thought in a real analysis
covariates_file_name = "Sample_Data/Covariates.txt"
output_file_name = "Sample_Data/eQTL_results_R.txt"