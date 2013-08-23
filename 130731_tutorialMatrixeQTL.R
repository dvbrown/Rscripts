library(MatrixEQTL)
setwd('~/Documents/eQTL/Matrix_eQTL_R/')
#This tutorial is hosted at http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/manual.html

#Matrix eQTL requires two input files, one with genotype information, and another with the gene expression. 
#It can also accept extra file with covariates and an error covariance matrix.

#The model of choice, there is more than 1 to choose from
useModel = modelLINEAR

#set up the SNP file and the expression file. The genotypes are coded as 0, 1 or 2 based on AA, AB, BB genotype.
#Column names must match
SNP_file_name = 'Sample_Data/SNP.txt'
expression_file_name = 'Sample_Data/GE.txt'

#The covariate input file. This will require much thought in a real analysis
covariates_file_name = "Sample_Data/Covariates.txt"
output_file_name = "Sample_Data/eQTL_results_R.txt"

#set this to 1e-5
pvOutputThreshold = 1e-2
#a dummy variable that is not often used
errorCovariance = numeric()

snps = SlicedData$new()
snps$fileDelimiter = "\t" # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values;
snps$fileSkipRows = 1 # one row of column labels
snps$fileSkipColumns = 1 # one column of row labels
snps$fileSliceSize = 2000 # read file in pieces of 2,000 rows
snps$LoadFile(SNP_file_name)

gene = SlicedData$new()
gene$fileDelimiter = "\t" # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1 # one row of column labels
gene$fileSkipColumns = 1 # one column of row labels
gene$fileSliceSize = 2000 # read file in pieces of 2,000 rows
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new()

#the main meat of the script. both snps and gene are the matrix objects that will be tested.
me = Matrix_eQTL_engine(snps,
                        gene,
                        cvrt,
                        output_file_name,
                        pvOutputThreshold,
                        useModel,
                        errorCovariance,
                        verbose = TRUE,
                        pvalue.hist = 'qqplot') #change this parameter to change the plot

plot(me)
output = read.delim(file='Sample_Data/eQTL_results_M.txt')

#next time try out the separate analysis of cis and trans eQTLs. Need location parameter files.
