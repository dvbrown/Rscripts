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


#### Plot some metric generated from 96 well plate into a 96 well matrix ####
# Insert well positions
wellMap = data.frame (rown = rep (letters[1:8], each=12), coln = rep (1:12, 8))
platelay = cbind(wellMap, dat)

# Change the aesthetics according to dataset
ggplot(platelay, aes(y = factor(rown, rev(levels(rown))),x = factor(coln))) + 
  geom_point(aes(colour = GenesDetected, size= LibSize))  + theme_bw() + scale_size(range = c(5, 20)) +
  labs(x=NULL, y = NULL) + scale_color_gradient(low="blue", high="red")
