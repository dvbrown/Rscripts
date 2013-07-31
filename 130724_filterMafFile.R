#A script to filter the maf file into a vcf file? Some file that is good for the tools.

mutations = read.delim('~/Documents/public-datasets/firehose/mafDashboard/2013_06_27/test.maf', skip=4)

mutFilter = mutations[,c(1,)]