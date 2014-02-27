
setwd('~/Documents/ELDA/140226_2ndyearReviewSummary/')
data = read.delim('140226_sphereEficiency.txt')
pvals = read.delim('140226_differenceTests.txt')

# Add a column with stars describing if a test is significant
cData$star <- " "
cData$star[cData$adj_pVal < .05]  = "*"
cData$star[cData$adj_pVal < .01]  <- "**"
cData$star[cData$adj_pVal < .001] <- "***"

ggplot(cData, aes(x=gene.x, y=mean, fill=cd133)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    xlab("Gene") +
    ylab("Expression normalised to GAPDH") +
    scale_fill_hue(name="CD133")+#, Legend label, use darker colors
    ggtitle("Differential expression relative to CD133 negative") +
    scale_y_continuous(breaks=0:20*4) +
    # Setting vjust to a negative number moves the asterix up a little bit to make the graph prettier
    geom_text(aes(label=star), colour="black", vjust=-2, size=10) +
    theme_bw(base_size=16)