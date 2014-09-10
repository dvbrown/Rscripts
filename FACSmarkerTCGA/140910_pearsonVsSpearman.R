getwd()
library(WGCNA)
library(sqldf)
library(ggplot2)
source('~/Documents/Rscripts/multiplot.R')

function (geneExpressionMatrix = dat, gene='PROM1') {
    corrPval = corAndPvalue(x=geneExpressionMatrix[,gene], y=geneExpressionMatrix)
    # Extract the correlation and p-value from the returned list
    correlation = corrPval$cor
    # Measure the  coefficient of determination
    coeffDeter = correlation^2
    pVal = corrPval$p
    fdr = p.adjust(pVal, method='fdr')
    
    result = t(rbind(correlation, coeffDeter, pVal, fdr))
    colnames(result) = c('correlation', ' Rsquared', 'p-value', 'FDR')
    return (result)
}

source('~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R')
agilent = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt", row.names=1)
rnaseq = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix", row.names=1)
setwd('/Users/d.brown6/Documents/public-datasets/cancerBrowser/deDupAgilent/results/')
db <- dbConnect(SQLite(), dbname='~/Documents/public-datasets/cancerBrowser/deDupAgilent/coexpression.sqlite')

cd133Ag = dbReadTable(db, "cd133Allgenes")
cd44Ag$genes = row.names(cd44Ag)
cd44Ag = dbReadTable(db, "cd44Allgenes")

resultSpear = cor(agilent[,'PROM1'], agilent, method = "spearman", use ="pairwise.complete.obs")
resultPear = cor(agilent[,'PROM1'], agilent, method='pearson', use ="pairwise.complete.obs")
scatterCD133 = as.data.frame(cbind(as.vector(resultSpear), as.vector(resultPear)))
colnames(scatterCD133) = c('spearman', 'pearson')

# Scatter plot the 2 
cd133 = ggplot(data=scatterCD133, aes(x=V1, y=V2)) + 
            geom_point(shape=19, alpha=1/4) + geom_smooth(method=lm, colour='red') +
            xlab("Spearman") + ylab("Pearson") +
            ggtitle("Comparison of Spearman and Pearson correlation for CD133") +  # Set title
            theme_bw(base_size=18)

############## Do the same scatter for CD44 ####################
resultSpear = cor(agilent[,'CD44'], agilent, use ="pairwise.complete.obs", method = "spearman")
resultPear = cor(agilent[,'CD44'], agilent, method='pearson', use ="pairwise.complete.obs")
scatterCD44 = as.data.frame(cbind(as.vector(resultSpear), as.vector(resultPear)))
colnames(scatterCD44) = c('spearman', 'pearson')

# Scatter plot the 2 
cd44 = ggplot(data=scatterCD44, aes(x=V1, y=V2)) + 
            geom_point(shape=19, alpha=1/4) + geom_smooth(method=lm, colour='red') +
            xlab("Spearman") + ylab("Pearson") +
            ggtitle("Comparison of Spearman and Pearson correlation for CD44") +  # Set title
            theme_bw(base_size=18)

multiplot(cd133, cd44, cols=1)
dbDisconnect(db) 

# Write out data
setwd('~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/spearmanPearsonComparison/correlations/')
database <- dbConnect(SQLite(), dbname='correlations.sqlite')
dbWriteTable(conn = database, name = "cd133Correlation", value = scatterCD133, row.names = TRUE)
dbWriteTable(conn = database, name = "cd44correlation", value = scatterCD44, row.names = TRUE)
dbDisconnect(database) 