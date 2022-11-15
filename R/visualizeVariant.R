#' Visualizes EVE score for a protein variant.
#'
#' A function that visualizes EVE scores of a variant at residue positions where
#' EVE scores are possible.  Function will replace EVE scores of NaN values with 0.
#' The colors assigned to EVE scores represent pathogenicity with 1 as most
#' pathogenic and 0 as benign.
#'
#' @param eveInfo Tibble with EVE scores for each residue position that has a
#' score calculated by EVE, residue position, wildtype amino acid, and mutated
#' amino acid.  If there are NaNs it means the variant provided doesn't have an
#' EVE score for that particular mutation at that position.  These NaN values
#' will be replaced by 0.
#'
#' @param geneName A character vector that is the name of the gene.  If not
#' supplied the default is X.
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
#' @importFrom methods hasArg
#' @import ggplot2 dplyr

visualizeVariant <- function(eveInfo, geneName = "X", aboveZeroOnly = FALSE) {

  # check if eveInfo has correct columns
  eveInfoCols <- colnames(eveInfo)
  expectedEveCols <- c("resPos", "eveScores", "wtAa", "varAa")

  if (!dplyr::setequal(eveInfoCols, expectedEveCols)) {
    stop("eveInfo does not have expected columns.  Ensure you used
         getEveScores to get the EVE scores for this variant first.")
  }

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

  p <- ggplot2::ggplot(eveInfoCopy, aes(x=resPos, y=eveScores)) +
    ggplot2::geom_segment(aes(x=resPos, xend=resPos, y=0, yend=eveScores),
                          color="grey") +
    ggplot2::geom_point(aes(color=eveScores), size=4) +
    ggplot2::scale_colour_gradient2(
      low = "steelblue1",
      mid = "gray",
      high = "firebrick1",
      midpoint = 0.5) +
    ggplot2::theme_light() +
    ggplot2::theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    ggplot2::labs(title = paste("EVE scores vs Residue Positions for", geneName, sep = " "),
         x = "Residue Position",
         y = "EVE Score",
         color = "EVE Score")

  return(p)
}

#' Visualize EVE scores of two variants for one gene simultaneously.
#'
#' A function that visualizes EVE scores for two variants of the same gene
#' simultaneously.  Function will replace EVE scores of NaN values with 0.
#'
#' @param eveInfo1 Tibble for the first variant with EVE scores for each residue
#' position that has a score calculated by EVE, residue position, wildtype amino
#' acid, and mutated amino acid.  If there are NaNs it means the variant
#' provided doesn't have an EVE score for that particular mutation at that
#' position.  These NaN values will be replaced by 0.
#'
#' @param eveInfo2 Tibble for the second variant with EVE scores for each
#' residue position that has a score calculated by EVE, residue position,
#' wildtype amino acid, and mutated amino acid.  If there are NaNs it means the
#' variant provided doesn't have an EVE score for that particular mutation at
#' that position.  These NaN values will be replaced by 0.
#'
#' @param geneName A character vector that is the name of the gene.  If not
#' supplied the default is X.
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
#' @importFrom methods hasArg
#' @import ggplot2 dplyr

visualizeVariant2 <- function(eveInfo1, eveInfo2, geneName = "X", aboveZeroOnly = FALSE) {
  # this one compares variants of the same gene, for example variants from
  # different samples

  # check if eveInfo1 and eveInfo2 have correct columns
  eveInfo1Cols <- colnames(eveInfo1)
  eveInfo2Cols <- colnames(eveInfo2)
  expectedEveCols <- c("resPos", "eveScores", "wtAa", "varAa")

  if (!dplyr::setequal(eveInfo1Cols, expectedEveCols)) {
    stop("eveInfo1 does not have expected columns.  Ensure you used
         getEveScores to get the EVE scores for this variant first.")
  }

  if (!dplyr::setequal(eveInfo2Cols, expectedEveCols)) {
    stop("eveInfo2 does not have expected columns.  Ensure you used
         getEveScores to get the EVE scores for this variant first.")
  }

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

  p <- ggplot2::ggplot(eveInfo1Copy, aes(x=resPos, y=eveScores)) +
    ggplot2::geom_segment(data=eveInfo1Copy, aes(x=resPos, xend=resPos, y=0, yend=eveScores),
                          color="grey") +
    ggplot2::geom_point(data=eveInfo1Copy, aes(color="Variant 1"), size=4) +
    ggplot2::geom_segment(data=eveInfo2Copy, aes(x=resPos, xend=resPos, y=0, yend=eveScores),
                          color="grey") +
    ggplot2::geom_point(data=eveInfo2Copy, aes(color="Variant 2"), size=4) +
    ggplot2::theme_light() +
    ggplot2::theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    ggplot2::labs(title = paste("EVE scores vs Residue Positions for", geneName, sep = " "),
         x = "Residue Position",
         y = "EVE Score",
         color = "EVE Score")

  return(p)
}

# [END]
