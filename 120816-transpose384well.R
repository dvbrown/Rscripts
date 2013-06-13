setwd('~/Documents/RNA data/qPCRexpt/')
map = read.delim('./120716-newPrimer/120816-newPrimer.txt', header=FALSE)
#transpose = vector(mode='character', length=384)

index = 1
while (index <= 9) {
  primers = map[index,]
  result = append(transpose, primers)
  index = index+1
}