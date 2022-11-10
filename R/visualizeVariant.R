#' Visualizes EVE score for a protein variant.
#'
#' A function that visualizes EVE scores of a variant at residue positions where
#' EVE scores are possible.
#'
#' @param eveInfo Tibble with EVE scores for each residue
#' position that has a score calculated by EVE, residue position, wildtype amino
#' acid, and mutated amino acid.  If there are NaNs it means the variant
#' provided doesn't have an EVE score.
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
#' variantPlotProt <- visualizeVariant(eveScoresProt)
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
#' variantPlotGen <- visualizeVariant(eveScoresGen)
#' variantPlotGen
#'
#' @import ggplot2

visualizeVariant <- function(eveInfo) {
  # TODO: add check to make sure has those columns
  # TODO: handle the NaN values and turn to 0
  p <- ggplot(eveInfo, aes(x=resPos, y=eveScores)) +
    geom_segment( aes(x=resPos, xend=resPos, y=0, yend=eveScores), color="grey") +
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
    labs(x = "Residue Position", y = "EVE Score", color = "EVE Score")
  return(p)
}

#' Visualize EVE scores of two variants for one gene.
#'
#' A function that visualizes EVE scores for two variants of the same gene
#' simultaneously.
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
#' compareVariantsPlot <- visualizeVariant2(eveScoresProt, eveScoresGen)
#' compareVariantsPlot
#'
#' @import ggplot2

visualizeVariant2 <- function(eveInfo1, eveInfo2) {
  # this one compares variants from different samples onto one gene
  p <- ggplot(eveInfo1, aes(x=resPos, y=eveScores)) +
    geom_segment(data=eveInfo1, aes(x=resPos, xend=resPos, y=0, yend=eveScores), color="grey") +
    geom_point(data=eveInfo1, aes(color="Variant 1"), size=4) +
    geom_segment(data=eveInfo2, aes(x=resPos, xend=resPos, y=0, yend=eveScores), color="grey") +
    geom_point(data=eveInfo2, aes(color="Variant 2"), size=4) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(x = "Residue Position", y = "EVE Score", color = "EVE Score")
  return(p)
}

# [END]
