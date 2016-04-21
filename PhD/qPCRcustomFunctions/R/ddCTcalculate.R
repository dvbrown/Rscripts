ddCTcalculate <-
function(geneOfInterest, sampleOfInterest='020_N', houseKeepingGene='GAPDH', referenceSample='020_N', data=rawData) {
  # The gene of interest and the sample of interest are both vectors from the dataFrame.
  # House keeping gene and 
  
  # This function has been successfully debugged and has been shown to work.
  sampleHouse = paste(sampleOfInterest, houseKeepingGene)
  SampleGene = paste(sampleOfInterest, geneOfInterest)
  # Extract the Cp of the hosue keeping gene
  houseCp = data[sampleHouse, 'meanCP']
  # place holder 020_N ATP5G3 as the current row
  geneCp = data[SampleGene, 'meanCP']
  # dCt calculation for the sample of interest
  dCt = houseCp - geneCp
  
  # Extract the meanCP for the reference sample. First get the index of the housekeeping gene, then the gene of interest
  refDctRowHouse = paste(referenceSample, houseKeepingGene)
  refDctRowGene = paste(referenceSample, geneOfInterest)
  # Calculate dCt for the reference sample
  referenceSample_dCt = data[refDctRowHouse, 'meanCP'] - data[refDctRowGene, 'meanCP']
  
  # Calculate ddCt
  ddCtNotSquared = dCt - referenceSample_dCt
  ddCt = 2^ddCtNotSquared
  
  return (ddCt)
}
