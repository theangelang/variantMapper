scoreVariant <- function(eveScores, weights = rep(1/length(eveScores), length(eveScores))) {
  # what to do if there are NaN, can replace with 0 and don't need to change the weights
  # or can change the weight so take that weight and redistribute

  # change weights name maybe since its a func too
  eveScoresCopy <- eveScores
  weightsCopy <- weights

  # check if eveScores and weights are the same length
  if (hasArg(weights)) {
    if (length(eveScores) != length(weights)) {
      warning("The length of the weights vector isn't the same as the length of eveScores.
              The default for weights will be used instead.")
      weightsCopy <- rep(1/length(eveScores), length(eveScores))
    }
  }

  # check if there are NaN values

  if(sum(is.nan(eveScores)) != 0) {
    warning("eveScores contains NaN values.  They will be replaced with 0 values.")
    eveScoresCopy[is.nan(eveScoresCopy)] <- 0
  }

  # calculate the weighted mean

  return(weighted.mean(eveScoresCopy, weightsCopy))

}
