library(ggplot2)
source('~/Code/Rscripts/Templates/multiplot.R')

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
correlation = ggplot(data=mungedData, aes(x=Cp.y, y=Cp.x, color=gene.x, label=row.names(mungedData))) + 
  geom_point(shape=19) + geom_smooth(method=lm, colour='red') + geom_text() +
  xlab("Replicate 1") + ylab("Replicate 2") + # Set axis labels
  ggtitle("Correlation of technical replicates") +  # Set title
  geom_hline(yintercept=13, colour = "green") + # Set a vertical line
  theme_bw(base_size=18)

#### Plot a histogram
p <- ggplot(df, aes(factor(category), dataPoints)) +
  geom_boxplot() + geom_jitter() +
  ggtitle("title") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
p

#### Plot a violin which is nicer than dotplot
p5 <- ggplot(dfLong_filter, aes(x=hashtag, y=log2(value+1), colour = factor(experiment))) +
  geom_violin() + geom_jitter(height = 0, width = 0.2) +
  ggtitle("Log hashtag counts") +
  xlab("hashtag") + ylab("log2 counts") +
  theme_bw(base_size=20) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p5

#### Make a dot plot vertical stacking of multiple groups
p1 <- ggplot(df, aes(factor(chr), counts, fill = factor(groups), label = factor(treatment))) +
    geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge") +
    scale_fill_manual(values=c("blue", "red"), name="Lysis buffer",    # Change the fill colour manually and set label title
                      breaks=c("group_1", "group_2"), labels=c("Group A", "GROUP B")) +     # Set the labels by first mapping to the groups then replacing them
    ggtitle("Percent mtDNA ATAC-seq reads") + geom_text(check_overlap = TRUE, size = 3) +
    xlab("% mtDNA reads") + ylab("group") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=16))
p1

#### Make a stacked barchart ####

plt5 <- ggplot(stats.pbmc) +
  geom_bar(aes(y=P, x= plate, fill=type), data=stats.pbmc, stat="identity") +
  ggtitle("") +
  xlab("Plate") + ylab("Percent reads") + 
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=22))
plt5

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
