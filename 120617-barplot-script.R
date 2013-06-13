housePlot=  ggplot(tbpMean, aes(x=sample, y=mean.Cp)) +
              geom_bar(position=position_dodge(width=1), fill='blue', colour='black') +
              geom_errorbar(aes(ymin=mean.Cp-sd.Cp, ymax=mean.Cp+sd.Cp), 
                width=.2,
                position=position_dodge(.9), colour='black') +
                  ylab('Mean Ct') + opts(axis.text.x = theme_text(angle=90, hjust=1.2, size=12), 
                                         title='Crossing point housekeeping gene GAPDH, TBP and PKM2')