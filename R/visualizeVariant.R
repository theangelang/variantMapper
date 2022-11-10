#' Visualizes EVE score for a protein variant.
#'
#' A function that visualizes EVE scores of a variant at residue positions where
#' EVE scores are possible.  Function will replace EVE scores of NaN values with 0.
#'
#' @param eveInfo Tibble with EVE scores for each residue
#' position that has a score calculated by EVE, residue position, wildtype amino
#' acid, and mutated amino acid.  If there are NaNs it means the variant
#' provided doesn't have an EVE score.
#'
#' @param geneName String that is the name of the gene.  If not supplied the
#' default is X.
#'
#' @param aboveZeroOnly Logical value that indicates whether to include EVE scores
#' of 0.  The default is FALSE so EVE scores of 0 will be included.  Note if
#' there are NaN values they will be excluded as well since they will be changed
#' to 0 first.
#'
#' @return A lollipop graph showing the EVE score at each residue position.
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
#' variantPlotProt <- visualizeVariant(eveScoresProt, aboveZeroOnly = TRUE)
#' variantPlotProt
#'
#' # If the data is in genomic form.
#' varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
#' varDataGen <- processVariantData(varDataGenPath, protein = FALSE)
#' varDataGen
#'
#' eveScoresGen <- getEveScores(EveData, varDataGen, protein = FALSE)
#' eveScoresGen
#'
#' variantPlotGen <- visualizeVariant(eveScoresGen, "NRXN1")
#' variantPlotGen
#'
#' @import ggplot2

visualizeVariant <- function(eveInfo, geneName = "X", aboveZeroOnly = FALSE) {
  # TODO: add check to make sure has those columns

  eveInfoCopy <- eveInfo

  # Replace NaN with 0
  if(sum(is.nan(eveInfo$eveScores)) != 0) {
    warning("eveInfo contains NaN values. They will be replaced with 0 values.")
    eveInfoCopy$eveScores[is.nan(eveInfoCopy$eveScores)] <- 0
  }

  if (methods::hasArg(aboveZeroOnly)) {
    if (isTRUE(aboveZeroOnly)) {
      eveInfoCopy <- dplyr::filter(eveInfoCopy, eveScores != 0)
      }
    }

  p <- ggplot(eveInfoCopy, aes(x=resPos, y=eveScores)) +
    geom_segment(aes(x=resPos, xend=resPos, y=0, yend=eveScores), color="grey") +
    geom_point(aes(color=eveScores), size=4) +
    scale_colour_gradient2(
      low = "steelblue1",
      mid = "gray",
      high = "firebrick1",
      midpoint = 0.5) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(title = paste("EVE scores vs Residue Positions for", geneName, sep = " "),
         x = "Residue Position",
         y = "EVE Score",
         color = "EVE Score")

  return(p)
}

#' Visualize EVE scores of two variants for one gene.
#'
#' A function that visualizes EVE scores for two variants of the same gene
#' simultaneously.  Function will replace EVE scores of NaN values with 0
#'
#' @param eveInfo1 Tibble for the first variant with EVE scores for each residue
#' position that has a score calculated by EVE, residue position, wildtype amino
#' acid, and mutated amino acid.  If there are NaNs it means the variant
#' provided doesn't have an EVE score.
#'
#' @param eveInfo2 Tibble for the second variant with EVE scores for each
#' residue position that has a score calculated by EVE, residue position,
#' wildtype amino acid, and mutated amino acid.  If there are NaNs it means the
#' variant provided doesn't have an EVE score.
#'
#' @param geneName String that is the name of the gene.  If not supplied the
#' default is X.
#'
#' @param aboveZeroOnly Logical value that indicates whether to include EVE scores
#' of 0.  The default is FALSE so EVE scores of 0 will be included.  Note if
#' there are NaN values they will be excluded as well since they will be changed
#' to 0 first.
#'
#' @return A lollipop graph showing the EVE score at each residue position for
#' both variants.
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
#' # If the data is in genomic form.
#' varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
#' varDataGen <- processVariantData(varDataGenPath, protein = FALSE)
#' varDataGen
#'
#' eveScoresGen <- getEveScores(EveData, varDataGen, protein = FALSE)
#' eveScoresGen
#'
#' compareVariantsPlot <- visualizeVariant2(eveScoresProt, eveScoresGen, "NRXN1")
#' compareVariantsPlot
#'
#' compareVariantsPlotAboveZero <- visualizeVariant2(eveScoresProt,
#' eveScoresGen, "NRXN1", aboveZeroOnly = TRUE)
#' compareVariantsPlotAboveZero
#'
#' @import ggplot2

visualizeVariant2 <- function(eveInfo1, eveInfo2, geneName = "X", aboveZeroOnly = FALSE) {
  # this one compares variants from different samples onto one gene
  # TODO: add check to make sure has those columns

  eveInfo1Copy <- eveInfo1
  eveInfo2Copy <- eveInfo2

  # Replace NaN with 0
  if(sum(is.nan(eveInfo1$eveScores)) != 0) {
    warning("eveInfo1 contains NaN values. They will be replaced with 0 values.")
    eveInfo1Copy$eveScores[is.nan(eveInfo1Copy$eveScores)] <- 0
  }

  if(sum(is.nan(eveInfo2$eveScores)) != 0) {
    warning("eveInfo2 contains NaN values. They will be replaced with 0 values.")
    eveInfo2Copy$eveScores[is.nan(eveInfo2Copy$eveScores)] <- 0
  }

  if (methods::hasArg(aboveZeroOnly)) {
    if (isTRUE(aboveZeroOnly)) {
      eveInfo1Copy <- dplyr::filter(eveInfo1Copy, eveScores != 0)
      eveInfo2Copy <- dplyr::filter(eveInfo2Copy, eveScores != 0)
      }
    }

  p <- ggplot(eveInfo1Copy, aes(x=resPos, y=eveScores)) +
    geom_segment(data=eveInfo1Copy, aes(x=resPos, xend=resPos, y=0, yend=eveScores), color="grey") +
    geom_point(data=eveInfo1Copy, aes(color="Variant 1"), size=4) +
    geom_segment(data=eveInfo2Copy, aes(x=resPos, xend=resPos, y=0, yend=eveScores), color="grey") +
    geom_point(data=eveInfo2Copy, aes(color="Variant 2"), size=4) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(title = paste("EVE scores vs Residue Positions for", geneName, sep = " "),
         x = "Residue Position",
         y = "EVE Score",
         color = "EVE Score")

  return(p)
}

# [END]
