getwd()
library(WGCNA)
library(ggplot2)
source('/Users/d.brown6/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R')
source('/Users/d.brown6/Documents/Rscripts/multiplot.R')
list.files()

dat = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt", row.names=1)
setwd('/Users/d.brown6/Documents/public-datasets/cancerBrowser/deDupAgilent/results/')

######################################## CD133 coexpressed Genes ################################################
cd133 = correlateGeneWithGEM(dat, 'PROM1')
# write.table(cd133, './140526_cd133Coexpression.txt', row.names=T, sep='\t')
plotCoexpression(cd133, 'CD133')

# Use twice the standard deviation and significantly correlated
cd133genes = cd133[(cd133[,1]) > 2*sd(cd133[,1]) & cd133[,4] < 0.05,]
# write.table(cd133genes, './140527_cd133Cutoff.txt', sep='\t')
cd133Square = makeSquareCoexpressionMatrix(cd133genes, dat)

cd133Dissim = makeDissimilarity(cd133Square)

######################################## CD44 coexpressed Genes ################################################
cd44 = correlateGeneWithGEM(dat, 'CD44')
plotCoexpression(cd44, 'CD44')

# Subset the dataframe with correlation values for those with high correlation and significance
# Use twice the standard deviation and significantly correlated
cd44genes = cd44[cd44[,1] > 2*sd(cd44[,1]) & cd44[,4] < 0.05,]

cd44Square = makeSquareCoexpressionMatrix(cd44genes, dat)
cd44Dissim = makeDissimilarity(cd44Square)

######################################## CD15 coexpressed Genes ################################################
setwd('/Users/d.brown6/Documents/public-datasets/cancerBrowser/deDupAgilent/results/CD15')
cd15 = correlateGeneWithGEM(dat, 'FUT4')

plotCoexpression(cd15, 'CD15')

# Subset the dataframe with correlation values for those with high correlation and significance
# Use twice the standard deviation and significantly correlated
cd15genes = cd15[cd15[,1] > 2*sd(cd15[,1]) & cd15[,4] < 0.05,]

cd15Square = makeSquareCoexpressionMatrix(cd15genes, dat)
cd15Dissim = makeDissimilarity(cd15Square)

######################################## Bind correlation and FDR ################################################
correlations = cbind(cd133[,1], cd44[,1], cd15[,1])
fdrs = cbind(cd133[,4], cd44[,4], cd15[,4])
colnames(correlations) = c('CD133', 'CD44', 'CD15')
colnames(fdrs) = c('CD133', 'CD44', 'CD15')
correlations = as.data.frame(correlations)

######################################## Plot the doubles ################################################
# get the genes which have been thresholded
cd133Corr = row.names(cd133genes)
cd44Corr = row.names(cd44genes)
cd15Corr = row.names(cd15genes)
cd44Cd15 = intersect(cd44Corr, cd15Corr)

# Set the colours of the scatterplots
correlations$threshold = "not significant"
correlations[cd133Corr,]$threshold = "CD133"
correlations[cd44Corr,]$threshold = "CD44"
correlations[cd15Corr,]$threshold = "CD15"
correlations[cd44Cd15,]$threshold = "CD44 + CD15"


######################################## Plot the doubles ################################################

cd133_44 = ggplot(data=correlations, aes(x=CD133, y=CD44, color=threshold)) + 
            geom_point(shape=19, alpha=1/4) + geom_smooth(method=lm, colour='red') +
            scale_colour_manual(values=c('red','green', 'blue', 'purple','grey')) +    
            xlab("CD133") + ylab("CD44") + # Set axis labels
            ggtitle("Correlation of genes with double FACS markers") +  # Set title
            theme_bw(base_size=18)

cd133_15 = ggplot(data=correlations, aes(x=CD133, y=CD15, color=threshold)) + 
    geom_point(shape=19, alpha=1/4) + geom_smooth(method=lm, colour='red') +
    scale_colour_manual(values=c('red','green', 'blue', 'purple','grey')) +
    xlab("CD133") + ylab("CD15") + # Set axis labels
    ggtitle("Correlation of genes with double FACS markers") +  # Set title
    theme_bw(base_size=18)

cd44_15 = ggplot(data=correlations, aes(x=CD44, y=CD15, color=threshold)) + 
    geom_point(shape=19, alpha=1/4) + geom_smooth(method=lm, colour='red') +
    scale_colour_manual(values=c('red','green', 'blue', 'purple','grey')) +
    xlab("CD44") + ylab("CD15") + # Set axis labels
    ggtitle("Correlation of genes with double FACS markers") +  # Set title
    theme_bw(base_size=18)

multiplot(cd133_15, cd133_44, cd44_15, cols=2)