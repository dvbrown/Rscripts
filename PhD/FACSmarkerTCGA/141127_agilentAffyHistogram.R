# A script to demonstrate Agilent is better than Affymetrix
library(ggplot2)
source("~/Documents/Rscripts/multiplot.R")

#Load Agilent
setwd('~/Documents/public-datasets/cancerBrowser/TCGA_GBM_G4502A_07_2-2014-05-02/')
agilent = read.delim("genomicMatrix", row.names=1)

# Load Affy
affy = read.delim('~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_u133a-2014-08-22/genomicMatrix', row.names=1)

patients = intersect(colnames(affy), colnames(agilent))
# Take patients in common
agilent = agilent[,patients]
affy = affy[,patients]

affyGAPDH = as.data.frame(t(affy["GAPDH",]))
affyGAPDH$gene = "GAPDH"
affyBact = as.data.frame(t(affy["ACTB",]))
affyBact$gene = "ACTB"
colnames(affyGAPDH) = c("expression", "gene")
colnames(affyBact) = c("expression", "gene")
affyGene = rbind(affyGAPDH, affyBact)

agilentGAPDH = as.data.frame(t(agilent["GAPDH",]))
agilentGAPDH$gene = "GAPDH"
agilentBact = as.data.frame(t(agilent["ACTB",]))
agilentBact$gene = "ACTB"
colnames(agilentGAPDH) = c("expression", "gene")
colnames(agilentBact) = c("expression", "gene")
agilentGene = rbind(agilentGAPDH, agilentBact)

affyHist = ggplot(affyGene, aes(x=expression, fill=gene)) + geom_density(alpha=.3) +
    xlab("Expression") + ylab("Density") +
    ggtitle("Affymetrix") +  # Set title
    theme_bw(base_size=24)

agilenist = ggplot(agilentGene, aes(x=expression, fill=gene)) + geom_density(alpha=.3) +
            xlab("Expression") + ylab("Density") +
            ggtitle("Agilent") +  # Set title
            theme_bw(base_size=24)

affyQQ = ggplot(affyGene, aes(sample = expression,  colour = gene)) + geom_point(stat = "qq") +
            xlab("Expression") + ylab("Density") +
            ggtitle("Affymetrix") +  # Set title
            theme_bw(base_size=24)

agilentQQ = ggplot(agilentGene, aes(sample = expression,  colour = gene)) + geom_point(stat = "qq") +
                xlab("Expression") + ylab("Density") +
                ggtitle("Agilent") +  # Set title
                theme_bw(base_size=24)

multiplot(affyQQ, agilentQQ, cols=1)
multiplot(affyHist, agilenist, cols=1)
multiplot(affyQQ, agilentQQ, affyHist, agilenist, cols=2)

rm(affy, agilent)
