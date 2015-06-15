library(ggplot2)
library(reshape)
library(plyr)
setwd("150612_hyperHypoxia/")

subtractBaseline = function(dataFrame, baseline, value) {
    dataSorted = dataFrame[order(dataFrame[,'PDGC'], dataFrame[,'Treatment']),]
    baseline = dataSorted[dataSorted$Treatment %in% baseline,]
    treatment = dataSorted[dataSorted$Treatment %in% value,]
    subtract =  treatment[,c(3:6)] - baseline[,c(3:6)]
    result = cbind(treatment[,c(1,2)], subtract)
    return (result)
}

dat_all = read.csv("150612_hyperHypoxia.csv",row.names=1)
colnames(dat_all) = c("PDGC", "Treatment", "Viable", "CD44+/ CD133-",
                  "CD44+/ CD133+", "CD44-/ CD133+", "CD44-/ CD133-")

dat = dat_all[,c(1,2,4:7)]
dat = dat[!dat$Treatment %in% "isotype",]

# Calculate percent change
change = rbind(subtractBaseline(dat, "control", "h2o2"),
               subtractBaseline(dat, "control", "hypoxia"))

# Convert from wide to long
mDat = melt(change, id.vars = c("PDGC", "Treatment"))

# Biological replicates
bioRep = ddply(mDat, .(PDGC, Treatment, variable), summarise, meanValue = mean(value, na.rm=T), 
               sdDiff = sd(value, na.rm=T), reps=length(value))
bioRep$seDiff = bioRep$sdDiff / (bioRep$reps)
bioRep$number = rep(c(2,4,3,1), 4)

write.csv(bioRep, "150614_replicates.csv")

#### Barchart ####
# Plot the biological replicates
mu020 = ggplot(data=bioRep[bioRep$PDGC %in% "MU020",], 
              aes(x=number, y=meanValue, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    ggtitle("Biological Replicates n = 3") +  ylim(-20, 30) +
    scale_fill_manual(values=c("orange", "blue")) + 
    geom_errorbar(aes(ymin=meanValue-seDiff, ymax=meanValue+seDiff), width=.2, position=position_dodge(0.9)) +
    xlab("Subpopulation") + ylab("Percent difference \nrelative to control treatment") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
mu020

mu039 = ggplot(data=bioRep[bioRep$PDGC %in% "MU039",], 
               aes(x=number, y=meanValue, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    # scale_y_continuous(breaks = round(seq(-20, 20, by = 5),1)) +
    ylim(-20, 30) +
    ggtitle("Biological Replicates n = 3") +  
    scale_fill_manual(values=c("orange", "blue")) + 
    geom_errorbar(aes(ymin=meanValue-seDiff, ymax=meanValue+seDiff), width=.2, position=position_dodge(0.9)) +
    xlab("Subpopulation") + ylab("Percent difference \nrelative to control treatment") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
mu039 + theme(text = element_text(size=20))

### Boxplot ###
box = ggplot(data=mDat, aes(x=number, y=value)) +
    geom_boxplot() + geom_point(aes(colour=Treatment), size =3, alpha=0.7,  position = position_jitter(w = 0.175)) +
    ggtitle("qPCR Summary") +
    xlab("Gene") + ylab("ddCt relative to CD133 negative") +
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    theme(text = element_text(size=20))
box

#### Stats ####
testDf = melt(dat, id.vars = c("PDGC", "Treatment"))
test20 = dat[dat$PDGC %in% "MU020",]
t.test(`CD44+/ CD133-` ~ Treatment, data=test20 ,subset = !test20$Treatment == "h2o2") # 5.132e-05
t.test(`CD44+/ CD133+` ~ Treatment, data=test20 ,subset = !test20$Treatment == "h2o2") # 0.0009871
t.test(`CD44-/ CD133+` ~ Treatment, data=test20 ,subset = !test20$Treatment == "h2o2") # 9.372e-05
t.test(`CD44-/ CD133-` ~ Treatment, data=test20 ,subset = !test20$Treatment == "h2o2") #  0.000317

t.test(`CD44+/ CD133-` ~ Treatment, data=test20 ,subset = !test20$Treatment == "hypoxia") # 0.9
t.test(`CD44+/ CD133+` ~ Treatment, data=test20 ,subset = !test20$Treatment == "hypoxia") # 0.15
t.test(`CD44-/ CD133+` ~ Treatment, data=test20 ,subset = !test20$Treatment == "hypoxia") # 0.27
t.test(`CD44-/ CD133-` ~ Treatment, data=test20 ,subset = !test20$Treatment == "hypoxia") #  0.1055


test39 = dat[dat$PDGC %in% "MU039",]
t.test(`CD44+/ CD133-` ~ Treatment, data=test39 ,subset = !test39$Treatment == "h2o2") # 2.464e-06
t.test(`CD44+/ CD133+` ~ Treatment, data=test39 ,subset = !test39$Treatment == "h2o2") # 0.0001852
t.test(`CD44-/ CD133+` ~ Treatment, data=test39 ,subset = !test39$Treatment == "h2o2") # 0.002159
t.test(`CD44-/ CD133-` ~ Treatment, data=test39 ,subset = !test39$Treatment == "h2o2") # 0.0007407

t.test(`CD44+/ CD133-` ~ Treatment, data=test39 ,subset = !test39$Treatment == "hypoxia") # 0.09
t.test(`CD44+/ CD133+` ~ Treatment, data=test39 ,subset = !test39$Treatment == "hypoxia") # 0.026
t.test(`CD44-/ CD133+` ~ Treatment, data=test39 ,subset = !test39$Treatment == "hypoxia") # 0.497
t.test(`CD44-/ CD133-` ~ Treatment, data=test39 ,subset = !test39$Treatment == "hypoxia") # 0.070

tet = read.csv("150614_replicatesPval.csv", row.names=1)
# Separate MU020 and MU039
tet20 = tet[tet$PDGC %in% "MU020",]
tet39 = tet[tet$PDGC %in% "MU039",]

tet20$correct = p.adjust(tet20$P.val, method = 'bonferroni')
tet39$correct = p.adjust(tet39$P.val, method = 'bonferroni')
tet = rbind(tet20, tet39)

write.csv(tet, "150614_replicatesPval.csv")