buildClassifier <-
function(signatureSurvivalFrame, percentileNum) {
    # Take a dataframe containing a signature score and censorship status and add the group membership eg 'high' or 'low'
    # Allows one to vary the percentile used as the classifier
    percent = quantile(signatureSurvivalFrame$sigScore, probs=percentileNum, names=T)
    signatureSurvivalFrame$percentile = ifelse(signatureSurvivalFrame$sigScore >= percent, 'high', 'low')
    return (signatureSurvivalFrame)
}
