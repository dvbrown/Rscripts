#correct the Agilent processed data in Natsume

#get the foreground and background intensities
data = read.delim('processed+data_expression_array-formatted2.txt', header=T)
background = subset(data, X == 'DarkCorner')
saturate = subset(data, X == 'GE_BrightCorner')
data.strip1 = subset(data, X != 'DarkCorner')
data.strip2 = subset(data.strip1, X != 'GE_BrightCorner')
back.mean = colMeans(background)
fore.mean = colMeans(saturate)

#Perform background correction
back.correct = apply(data.num, 2, "-", back.mean)
fore.correct = apply(back.correct, 2, "/", fore.mean )

#log transform data
log.data = apply(merge.data[2:14], 2, 'log2')