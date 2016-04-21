niceErrorBarPlot <-
function (summarisedData, xAxis=gene.x, yAxis=mean, groupVariable=cd133, 
                              title='A title', xLabel='Gene', yLabel='Expression') {
  p = ggplot(summarisedData, aes(x=xAxis, y=yAxis, fill=groupVariable)) + 
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=yAxis-se, ymax=yAxis+se),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9)) +
        xlab(xLabel) +
        ylab(yLabel) +
        scale_fill_hue(name="CD133")+#, Legend label, use darker colors
        ggtitle(title) +
        scale_y_continuous(breaks=0:20*4) +
        # Setting vjust to a negative number moves the asterix up a little bit to make the graph prettier
        geom_text(aes(label=star), colour="black", vjust=-2, size=10) +
        theme_bw(base_size=16)
  return (p)
}
