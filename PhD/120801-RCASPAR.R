#Match up the paitent numbers in the expression matrix to the survival matrix
sur.sub = surData[,c(4,5)]
expr.names = row.names(t.frame)
sur.sub = sur.sub[expr.names,]
#Get rid of the last paitent that only has NAs in expression
t.2 = t.frame[1:91,]
sur.2 = sur.sub[1:91,]

colnames(surData) = c('age.diagnosis','gender', 'karnofsky', 'True_STs', 'censored')
model = STpredictor_xvBLH(t.2, sur.2, k=5, cut.off=400)
model
#cannot get even partition into long and short term survivors