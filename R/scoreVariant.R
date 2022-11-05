scoreVariant <- function(eveScores, posWeights = rep(1/length(eveScores), length(eveScores))) {
  # what to do if there are NaN, can replace with 0 and don't need to change the weights
  # or can change the weight so take that weight and redistribute

  eveScoresCopy <- eveScores
  posWeightsCopy <- posWeights

  # check if eveScores and posWeights are the same length
  if (hasArg(posWeights)) {
    if (length(eveScores) != length(posWeights)) {
      warning("The length of the posWeights vector isn't the same as the length of eveScores.
              The default for posWeights will be used instead.")
      posWeightsCopy <- rep(1/length(eveScores), length(eveScores))
    }
  }

  # check if there are NaN values

  if(sum(is.nan(eveScores)) != 0) {
    warning("eveScores contains NaN values.  They will be replaced with 0 values.")
    eveScoresCopy[is.nan(eveScoresCopy)] <- 0
  }

  # calculate the weighted mean

  return(weighted.mean(eveScoresCopy, posWeightsCopy))

}
