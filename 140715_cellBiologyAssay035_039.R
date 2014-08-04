source('~/Documents/Rscripts/cellBiologyAnalysisFunctions.R')

############################################## IO ###############################################
setwd('~/Documents/Cell_biology/proliferation/Resazurin/140710_039_035/')
growthD3 = read.delim('140709_day3Rep.txt')
growthD7 = read.delim('140714_day7Rep.txt')
growthD3$clone = as.factor(growthD3$clone)
growthD7$clone = as.factor(growthD7$clone)
growthD3$subPop = as.factor(growthD3$subPop)
growthD7$subPop = as.factor(growthD7$subPop)

############################################## means and averages ###############################################

# Subtract background and take mean and SD
day3Growth = backgroundMeanSD(growthD3)
day7Growth = backgroundMeanSD(growthD7)

# Did not plate down the doublr positive
day3Growth = day3Growth[!day3Growth$subPop %in% 'CD44+/CD133+',]
day7Growth = day7Growth[!day7Growth$subPop %in% 'CD44+/CD133+',]

#### Plot the raw results ####
growthPlot3 = ggplot(data=day3Growth[day3Growth$treatment %in% 'DMSO',], 
                     aes(x=clone, y=mean, fill=subPop)) + 
    scale_fill_manual(values=c("darkorange", "royalblue", "forestgreen", "red")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Fluorescent intensity") +
    ggtitle("Comparing growth at day 3 by \nmarker status") +  # Set title
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

growthPlot7 = ggplot(data=day7Growth[day7Growth$treatment %in% 'DMSO',], 
                     aes(x=clone, y=mean, fill=subPop)) + 
    scale_fill_manual(values=c("darkorange", "royalblue", "forestgreen")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Fluorescent intensity") +
    ggtitle("Comparing growth at day 7 \nby marker status") +  # Set title
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#multiplot(growthPlot3, growthPlot7)

############################################## Calculate and plot the DMSO corrected values #################################################
day3TMZ = calcDMSOcontrol(day3Growth)
day7TMZ = calcDMSOcontrol(day7Growth)
# 
tmzPlot3 = ggplot(data=day3TMZ, aes(x=clone, y=dmsoCorrected, fill=subPop)) + 
    scale_fill_manual(values=c("gold", "chartreuse4", "skyblue2")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("Clone") + ylab("Cell number relative to \nDMSO control") +
    ggtitle("Temozolomide sensitivty at day 3 by \nmarker status") +  # Set title
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

tmzPlot7 = ggplot(data=day7TMZ, aes(x=clone, y=dmsoCorrected, fill=subPop)) + 
    scale_fill_manual(values=c("gold", "chartreuse4", "skyblue2")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("Clone") + ylab("Cell number relative to \nDMSO control") +
    ggtitle("Temozolomide sensitivty at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#multiplot(tmzPlot3, tmzPlot7)

#multiplot(growthPlot3, growthPlot7, tmzPlot3, tmzPlot7, cols=2)
# write.table(day3TMZ, '140715_day3TMZprocessed.txt', sep='\t')
# write.table(day7TMZ, '140715_day7TMZprocessed.txt', sep='\t')