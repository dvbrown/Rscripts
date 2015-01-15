# Create a supplementary table for Coexpression study with big sigs and verhhak intersect

cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)
cd15 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/CD15/140528_cd15Cutoff.txt", row.names=1)
aldh1 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ALDH1/140528_ALDH1Cutoff.txt", row.names=1)
itag6 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ITGA6//140528_ITGA6Cutoff.txt", row.names=1)
l1cam = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/L1CAM/140528_L1CAMCutoff.txt", row.names=1)

# Read in the verhaak signatures
verhaakSig = read.delim('~/Documents/public-datasets/TCGA/classficationSignature/131022_danFixedTCGAsignature.txt')
verSigs = list(verhaakSig$Proneural, verhaakSig$Neural, verhaakSig$Classical, verhaakSig$Mesenchymal)
names(verSigs) = colnames(verhaakSig)

bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), "CD15" = row.names(cd15),
               "ALDH1"=row.names(aldh1), "ITGA6"=row.names(itag6), "L1CAM"=row.names(l1cam))
suppTable = as.data.frame(bigSigs)

rm(cd133Sig, cd44Sig, cd15, aldh1, itag6, l1cam)

# Intersect CD133 - PN
cd133PN = intersect(bigSigs$CD133, verhaakSig$Proneural)

# Intersect CD44 - MES
cd44MES = intersect(bigSigs$CD44, verhaakSig$Mesenchymal)

# Intersect CD15 - MES
cd15MES = intersect(bigSigs$CD15, verhaakSig$Mesenchymal)

write.table(cd133PN, "~/Documents/manuscripts/140610_inSilicoGBM/corrections/cd133PN.txt", sep="\t")
write.table(cd44MES, "~/Documents/manuscripts/140610_inSilicoGBM/corrections/cd44MES.txt", sep="\t")
write.table(cd15MES, "~/Documents/manuscripts/140610_inSilicoGBM/corrections/cd15MES.txt", sep="\t")