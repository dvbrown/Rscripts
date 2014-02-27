
setwd('~/Documents/ELDA/140226_2ndyearReviewSummary/')
data = read.delim('140226_sphereEficiency.txt')
pvals = read.delim('140226_differenceTests.txt')
pvals = pvals[c(1:6),]

pvals$adj_pVal = p.adjust(pvals$pVal, 'fdr')

x = ggplot(data, aes(x=Group, y=Estimate, fill=cd133)) + 
    scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    xlab("Patient clone") +
    ylab("Sphere forming efficiency (inverse)") +
    ggtitle("Repopulation efficency of CD133 sorted cells") +
    scale_y_continuous(breaks=0:20*4) +
    # Setting vjust to a negative number moves the asterix up a little bit to make the graph prettier
    theme_bw(base_size=16)

# rotate x axis label
x + theme(axis.text.x = element_text(angle = 90, hjust = 1))
