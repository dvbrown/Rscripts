# Compare veerhak and my FACS signatures
setwd('~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/signatureComparison/')

########################## Read in the signatures ################################
rnaseqGem = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix", row.names=1)

cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)
cd15 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/CD15/140528_cd15Cutoff.txt", row.names=1)
aldh1 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ALDH1/140528_ALDH1Cutoff.txt", row.names=1)
itag6 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ITGA6//140528_ITGA6Cutoff.txt", row.names=1)
l1cam = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/L1CAM/140528_L1CAMCutoff.txt", row.names=1)

fascSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), "CD15" = row.names(cd15),
               "ALDH1"=row.names(aldh1), "ITGA6"=row.names(itag6), "L1CAM"=row.names(l1cam))

tcgaSigs = read.delim('~/Documents/public-datasets/TCGA/classficationSignature/131022_danFixedTCGAsignature.txt')

########################## CD133 and proneural ########################## 
proCD133 = intersect(tcgaSigs$Proneural, row.names(cd133Sig))

# Use phyper the hypergeometric distribution to compute overlap
# phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
# q = size of overlap-1; m=number of upregulated genes in experiment #1; n=(total number of genes on platform-m); k=number of upregulated genes in experiment #2.

phyper(length(proCD133) - 1, length(row.names(cd133Sig)), 20000 - length(row.names(cd133Sig)), length(tcgaSigs$Proneural), lower.tail=FALSE, log.p=FALSE)

########################## CD44 and Mesenchymal ########################## 
mesCD44 = intersect(tcgaSigs$Mesenchymal, row.names(cd44Sig))

phyper(length(mesCD44) - 1, length(row.names(cd44Sig)), 20000 - length(row.names(cd44Sig)), length(tcgaSigs$Mesenchymal), lower.tail=FALSE, log.p=FALSE)
