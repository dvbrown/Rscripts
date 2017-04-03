library(ggplot2)

#### Barchart ####
plt = ggplot(data=bioRep[bioRep$PDGC %in% "MU020",], 
              aes(x=number, y=meanValue, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    ggtitle("Biological Replicates n = 3") +  ylim(-20, 30) +
    scale_fill_manual(values=c("orange", "blue")) + 
    geom_errorbar(aes(ymin=meanValue-seDiff, ymax=meanValue+seDiff), width=.2, position=position_dodge(0.9)) +
    xlab("Subpopulation") + ylab("Percent difference \nrelative to control treatment") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
plt

#### Plot correlation of the replicate ####
correlation = ggplot(data=mungedData, aes(x=Cp.y, y=Cp.x, color=gene.x)) + 
  geom_point(shape=19) + geom_smooth(method=lm, colour='red') +
  xlab("Replicate 1") + ylab("Replicate 2") + # Set axis labels
  ggtitle("Correlation of technical replicates") +  # Set title
  theme_bw(base_size=18)

#### Plot a histogram
p <- ggplot(df, aes(factor(category), dataPoints)) +
  geom_boxplot() + geom_jitter() +
  ggtitle("title") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
p


#### Plot some metric generated from 96 well plate into a 96 well matrix ####

plot96wellMap <- function (sampleSubset, plt_Title, variable_2_colour, variable_2_size) {
  # DATAFRAME    sampleSubset = the 96 samples you want to plot. This must be in order with well A12 = sample 12 and well B1 = sample 13
  # STRING    plt_Title = The title of your plot
  # STRING     variable_2_colour = the measurement you want mapped to colour
  # STRING     variable_2_size = the measurement you want mapped to circle size
  wellMap = data.frame (rown = rep (LETTERS[1:8], each=12), coln = rep (1:12, 8))
  platelay = cbind(wellMap, sampleSubset)
  # ggplot using well coordinates generate dabove
  plt = ggplot(platelay, aes(y = factor(rown, rev(levels(rown))),x = factor(coln))) + 
    geom_point(aes_string(colour = variable_2_colour, size= variable_2_size)) + scale_size(range = c(5, 10)) +
    labs(x=NULL, y = NULL) + scale_color_gradient(low="blue", high="red") + ggtitle(plt_Title) +
    theme_bw() + theme(axis.text=element_text(size=16))
  return(plt)
}
