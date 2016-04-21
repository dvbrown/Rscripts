censorData <-
function(survivalData) {
    # Takes patient dataframe
    status = survivalData[,'patient.followups.followup.vitalstatus']
    status = gsub('alive', 0, status)
    status = gsub('dead', 1, status)
    status = gsub('deceased', 1, status)
    # Check if the survival analysis can handle NAs for survival
    result = survivalData[,c(1,4,5)]
    censored = cbind(result, status)
    return (censored)
}
