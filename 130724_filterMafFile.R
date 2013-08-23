#A script to filter the maf file into a vcf file? Some file that is good for the tools.

mutations = read.delim('~/Documents/public-datasets/firehose/mafDashboard/2013_06_27/test.maf', skip=4)

mutFilter = mutations[,c(1,9,16,26)]
#rm(mutations)

position = mutations[,c(1,5,6,7)]
position$Start_position = as.integer(position$Start_position)
position$End_position = as.integer(position$End_position)
position$Average_position = as.integer((position$Start_position + position$End_position)/2)
positionFinal = position[,c(1,2,5)]

mutFilter$Variant_Classification = as.character(mutFilter$Variant_Classification)
#remove those mutations classed as silent
mutClass = mutFilter[(mutFilter$Variant_Classification != 'Silent'),]
write.table(mutClass,'./130806_filteredMutationsTCGA.txt',sep='\t',row.names=FALSE)