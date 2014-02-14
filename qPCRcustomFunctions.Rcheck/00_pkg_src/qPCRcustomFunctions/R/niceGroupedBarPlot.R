niceGroupedBarPlot <-
function (dataFrame, ddCt='ddCt', sampleOrigin="origin.x", gene="gene.x", graphTitle="A pretty plot") {
  ggplot(data=dataFrame, aes(x=sampleOrigin, y=ddCt, fill=gene)) + 
      geom_bar(stat="identity", position=position_dodge(), colour="black") + 
      scale_fill_hue(name="Gene") +      # Set legend title
      xlab("Sample") + ylab("ddCt") + # Set axis labels
      ggtitle(graphTitle) +  # Set title
      theme_bw(base_size=18)
}
