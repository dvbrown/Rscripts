library(ggplot2)
setwd('~/Documents/ELDA/140226_2ndyearReviewSummary/')
data = read.delim('140227_spherePercent.txt')
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
    ylab("Sphere forming efficiency (percent)") +
    ggtitle("Repopulation potential of CD133 sorted cells") +
    scale_y_continuous(breaks=0:20*4) +
    # Setting vjust to a negative number moves the asterix up a little bit to make the graph prettier
    theme_bw(base_size=20)

# rotate x axis label
x + theme(axis.text.x = element_text(angle = 90, hjust = 1))

########################################## Subset the data for primary and recurrent #############################
data = data[data$Group %in% c('#030_N', '#030_P', '#030a_N', '#030a_P'),]
x = ggplot(data, aes(x=Group, y=Estimate, fill=cd133)) + 
    scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    xlab("Patient clone") +
    ylab("Sphere forming efficiency (percent)") +
    ggtitle("Repopulation potential of CD133 sorted cells") +
    scale_y_continuous(breaks=0:20*4) +
    # Setting vjust to a negative number moves the asterix up a little bit to make the graph prettier
    theme_bw(base_size=20)

# rotate x axis label
x + theme(axis.text.x = element_text(angle = 90, hjust = 1))

########################################## Now some normalised data #############################
# I calculated this is excel

negative = c(9.009426939,1.931042598)
positive = c(14.25994723,3.452871847)

summary = rbind(negative, positive)
colnames(summary) = c('mean', 'sem')
summary = as.data.frame(summary)
summary$group = row.names(summary)

y = ggplot(summary, aes(x=row.names(summary), y=mean, fill=group)) +
    scale_fill_manual(values=c("pink", "cyan")) +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=(mean-sem), ymax=(mean+sem)),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    xlab("CD133 status") +
    ylab("Normalised sphere forming efficiency (percent)") +
    ggtitle("Repopulation potential of CD133 sorted cells") +
    scale_y_continuous(breaks=0:20*4) +
    # Setting vjust to a negative number moves the asterix up a little bit to make the graph prettier
    theme_bw(base_size=20)

# rotate x axis label
y + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# t-test for group difference
t.test(Estimate~cd133, data=data, paired=T)