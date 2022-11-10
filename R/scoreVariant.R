#' Scores a variant of a protein based on weighted average.
#'
#' @param eveScores A vector containing EVE scores for each residue position
#' that has a possible EVE
#' score.
#'
#' @param posWeights The weights to give each residue position when calculating
#' the weighted average of EVE scores for the protein variant.  The default is
#' using equal weights for each residue position i.e. 1/(number of residues with
#' EVE scores).  This is used when an argument is not passed into the posWeights
#' parameter.
#'
#' @return Returns a double with the weighted EVE score for the protein variant.
#'
#' @export
#'
#' @examples
#' # Examples:
#' # First process the EVE data and variant data.
#' EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
#' EveData <- processEveData(EvePath)
#' EveData
#'
#' # If the data is in protein form.
#' varDataProtPath <- system.file("extdata", "variant_data_protein.csv", package = "variantMapper")
#' varDataProt <- processVariantData(varDataProtPath, protein = TRUE)
#' varDataProt
#'
#' eveScoresProt <- getEveScores(EveData, varDataProt, protein = TRUE)
#' eveScoresProt
#'
#' # Get the weighted average EVE score to see overall how pathogenic this variant is
#' avgEveScore <- scoreVariant(eveScoresProt$eveScores) #using the default of uniform weights
#' avgEveScore
#'
#' # If the data is in genomic form.
#' varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
#' varDataGen <- processVariantData(varDataGenPath, protein = FALSE)
#' varDataGen
#'
#' eveScoresGen <- getEveScores(EveData, varDataGen, protein = FALSE)
#' eveScoresGen
#'
#' # Get the weighted average EVE score to see overall how pathogenic this variant is
#' # Assign more weight to the first few residue positions
#' weightsUnequal <- rep(1/(2 * length(eveScoresGen$eveScores)), length(eveScoresGen$eveScores))
#' weightsUnequal[1] <- 0.3
#' weightsUnequal[2] <- 0.2
#' weightedAvgEveScore <- scoreVariant(eveScoresGen$eveScores, weightsUnequal)
#' weightedAvgEveScore
#'
#' @importFrom methods hasArg
#' @importFrom stats weighted.mean

scoreVariant <- function(eveScores, posWeights = rep(1/length(eveScores), length(eveScores))) {
  # what to do if there are NaN, can replace with 0 and don't need to change the weights
  # or can change the weight so take that weight and redistribute

  eveScoresCopy <- eveScores
  posWeightsCopy <- posWeights

  # check if eveScores and posWeights are the same length
  if (methods::hasArg(posWeights)) {
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

  return(stats::weighted.mean(eveScoresCopy, posWeightsCopy))

}

# [END]
