#' Helper function to get a list of unique wild type amino acids at each residue
#' position.
#'
#' @param eveData A tibble with a column "wtAa" and "resPos" with 1 letter amino
#' acid code and residue position respectively.
#'
#' @return A tibble with wild type amino acid and residue position.
#'
#' @noRd

uniqueWtAaPos <- function(eveData) {
  # Note- wtAa and resPos are columns in eveData
  # TODO: Include check to make sure they are in the eveData
  wtAaPos <- dplyr::distinct(eveData, wtAa, resPos)
  return(wtAaPos)
}

#' Helper function to find out which nucleotide position of a codon is mutated.
#'
#' @param coordinates A vector of doubles indicating the start position in the
#' genome for each residue.
#'
#' @param varCoord An integer representing where in the genome the mutation
#' occurred.
#'
#' @return A vector of length 2 with the nucleotide position (1, 2, or 3) in the
#'  first position and the genomic coordinate of the first nucleotide in the
#'  corresponding codon.  If the varCoord is not in the coordinates it will
#'  return a vector of length 2 with values of NaN for both values.
#'
#' @noRd

findVariantPosition <- function(coordinates, varCoord) {
  if (varCoord %in% coordinates) {
    return(c(ntPos = 1, genomicCoord = varCoord))
  } else if ((varCoord - 1) %in% coordinates) {
    return(c(ntPos = 2, genomicCoord = varCoord - 1))
  } else if ((varCoord - 2) %in% coordinates) {
    return(c(ntPos = 3, genomicCoord = varCoord - 2))
  } else {
    # the genomic coordiante is not scored
    return(c(ntPos = NaN, genomicCoord = NaN))
  }
}

#' Helper function to construct the alternate codon sequence corresponding to
#' the single nucleotide variant (SNV).
#'
#' @param refData A tibble with three columns POS for the genomic position, REF
#' for the wild type codon, and resPos for the residue position.
#'
#' @param ntPos Which position in the codon the SNV occurred the first, second or
#' third.
#'
#' @param genomicCoord The genomic coordinate of the first nucleotide in the
#' corresponding codon for the SNV.
#'
#' @param varNt The nucelotide variant.
#'
#' @return A tibble with reference codon, alternate codon, and residue position.
#'
#' @noRd

constructAltSeq <- function(refData, ntPos, genomicCoord, varNt) {
  # TODO: POS is column in refData, check if it is present
  # TODO: check if REF and resPos present too
  refInfo <- dplyr::filter(refData, POS == genomicCoord) # get the matching row
  refCodon <- refInfo[1,]$REF
  resPos <- refInfo[1,]$resPos

  # produce the alternative codon seq
  altCodon <- refInfo[1,]$REF # initially assign it ref codon seq
  stringr::str_sub(altCodon, ntPos, ntPos) <- varNt # change to alternate codon
                                                    # seq
  return(tibble(refCodon = refCodon, altCodon = altCodon, resPos = resPos))
}

#' Gets the EVE score for the amino acids in variant of interest.
#'
#' A function that takes in EVE and variant data then will get the EVE score for
#' each amino acid in the variant of interest.  If the data is supplied in
#' protein form it will take the average of the EVE scores for that amino acid
#' mutation since EVE data is available down to the single nucleotide variation
#' level and multiple codons code for the same amino acid.
#'
#' @param eveData A tibble containing the EVE data that has already been
#' processed by the processEveData function.
#'
#' @param variantData A tibble containing the variant data that has already
#' been processed by the processVariantData function.  Assumption is that there
#' will only be one variant at each residue/genomic position and the wildtype
#' and variant amino acid are distinct from one another.  For protein data it
#' will be missense variants only and genomic data single nucleotide variants
#' only.  Won't try to score mutations that result in wildtype and variant amino
#' acid being the same.  Assume variants are independent of each other.
#'
#' @param protein Specifies whether the variant data is in protein or genomic
#' form.  By default it is set to TRUE meaning the default is protein form.
#'
#' @return Returns a tibble containing the EVE scores for each residue
#' position that has a score calculated by EVE, residue position, wildtype amino
#' acid, and mutated amino acid.  If there are NaNs it means the variant
#' provided doesn't have an EVE score.
#'
#' @export
#'
#' @examples
#' # Example:
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
#' eveScores <- getEveScores(EveData, varDataProt, protein = TRUE)
#' eveScores
#'
#' # If the data is in genomic form.
#' varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
#' varDataGen <- processVariantData(varDataGenPath, protein = FALSE)
#' varDataGen
#'
#' eveScores <- getEveScores(EveData, varDataGen, protein = FALSE)
#' eveScores
#'
#' @importFrom stats setNames
#' @import dplyr tibble

getEveScores <- function(eveData, variantData, protein = TRUE) {

  # eveData is the tibble
  # return a vector with the EVE score for all the proteins with EVE score

  # get pairs of distinct wtAa and residue positions
  wtAaPos <- uniqueWtAaPos(eveData)
  numResidues <- nrow(wtAaPos)
  eveScores <- vector("numeric", numResidues) # vector to return
  eveScores <- stats::setNames(eveScores, wtAaPos$resPos)
  wtAas <- wtAaPos$wtAa
  wtAas <- stats::setNames(wtAas, wtAaPos$resPos)
  varAas <- wtAas
  varAas <- stats::setNames(varAas, wtAaPos$resPos)

  # those with no mutation i.e. wildtype have EVE score 0 (benign), and those
  # variants who aren't scored are assigned NaN

  if (isTRUE(protein)) {
    # no translation required

    # TODO: do a null check for variants
    variants <- dplyr::inner_join(eveData, variantData) # variants with EVE data

    eveCol <- which(colnames(variants) == "EVE") # get column with EVE score
    varAaCol <- which(colnames(variants) == "varAa")

    # loop through the variant data
    for (i in 1:nrow(variantData)) {
      mut <- variantData[i,]
      pos <- as.character(mut$resPos) # get position where to put the score

      # get a subset of all the rows in EVE data with that wtAa, resPos, varAa
      # TODO: check if it has all those columns
      varSubset <- dplyr::filter(variants,
                                 wtAa == mut$wtAa,
                                 resPos == mut$resPos,
                                 varAa == mut$varAa)

      # score the mutation
      if (nrow(varSubset) == 0) {
        # variant isn't scored
        eveScores[[pos]] <- NaN
        varAas[[pos]] <- mut$varAa# udpate the varAa
      } else if (nrow(varSubset) == 1) {
        eveScores[[pos]] <- varSubset[[1, eveCol]]
        varAas[[pos]] <- varSubset[[1, varAaCol]]# udpate the varAa
      } else {
        # assign EVE score as the average
        # TODO: look into if all variants mutate to same aa but different codon
        # have same EVE score
        eveScores[[pos]] <- mean(varSubset$EVE)
        varAas[[pos]] <- varSubset[[1, varAaCol]]# udpate the varAa
      }
    }
    result <- tibble::tibble(eveScores = unname(eveScores),
                             resPos = wtAaPos$resPos,
                             wtAa = unname(wtAas),
                             varAa = unname(varAas))
    return(result)
  } else {
    # it is in genomic coordinates

    # unique genomic coordinates, wt seq, and residue position in EVE data
    # TODO: check it has all those cols first
    uniqueGenomicCoords <- dplyr::distinct(eveData, POS, REF, resPos)

    # data processing checked if all variants are single nucleotide variants
    # (SNVs), will assume data consists of just SNVs

    # get the nucleotide position in the codon and corresponding start genomic
    # coordinate for the codon
    # TODO: check it has all those cols first
    ntPosAndGenomicCoord <- sapply(variantData$start,
                                   findVariantPosition,
                                   coordinates = uniqueGenomicCoords$POS)

    ntPosAndGenomicCoord <- t(ntPosAndGenomicCoord) # transpose data
    ntPosAndGenomicCoord <- dplyr::as_tibble(ntPosAndGenomicCoord)

    refCodons <- c("character", nrow(variantData))
    altCodons <- c("character", nrow(variantData))
    resPos <- c("numeric", nrow(variantData))

    #construct the alternative allele for each variant
    for (i in 1:nrow(variantData)) {
      # TODO: check it has all those columns first
      converted <- constructAltSeq(refData = uniqueGenomicCoords,
                                   ntPos = ntPosAndGenomicCoord[i,]$ntPos,
                                   genomicCoord = ntPosAndGenomicCoord[i,]$genomicCoord,
                                   varNt = variantData[i,]$ALT)
      refCodons[i] <- converted[1,]$refCodon
      altCodons[i] <- converted[1,]$altCodon
      resPos[i] <- converted[1,]$resPos
    }

    resPos <- as.numeric(resPos)

    # combine variant data with refCodons, altCodons, resPos, and genomic
    # coordinate of first nt in codon
    mutatedVariantData <- variantData %>%
      tibble::add_column(refCodons = refCodons,
                         altCodons = altCodons,
                         resPos = resPos,
                         genomicCoord = ntPosAndGenomicCoord$genomicCoord,
                         ntPos = ntPosAndGenomicCoord$ntPos)

    # get EVE data for variants
    # TODO: check if it has those columns first
    variants <- dplyr::inner_join(eveData,
                                  mutatedVariantData,
                                  by = c("POS" = "genomicCoord",
                                         "REF" = "refCodons",
                                         "ALT" = "altCodons"))

    # get EVE score for variants
    for (i in 1:nrow(variantData)) {
      mut <- mutatedVariantData[i,]
      pos <- as.character(mut$resPos) # get position where to put the score

      # find EVE data that matches variant
      # TODO: check if it has all those columns first
      varSubset <- dplyr::filter(variants,
                                 POS == mut$genomicCoord,
                                 resPos == mut$resPos,
                                 REF == mut$refCodons,
                                 ALT == mut$altCodons)

      if (nrow(varSubset) == 0) {
        # variant isn't scored
        eveScores[[pos]] <- NaN
        varAas[[pos]] <- mut$varAa # udpate the varAa
      } else if (nrow(varSubset) == 1) {
        eveScores[[pos]] <- varSubset[1,]$EVE
        varAas[[pos]] <- varSubset[1,]$varAa # udpate the varAa
      } else {
        # assign EVE score as the average
        # TODO: look into if all variants mutate to same aa but different codon
        # have same EVE score
        eveScores[[pos]] <- mean(varSubset$EVE)
        varAas[[pos]] <- varSubset[1,]$varAa # udpate the varAa
      }
    }
    return(tibble(eveScores = unname(eveScores),
                  resPos = wtAaPos$resPos,
                  wtAa = unname(wtAas),
                  varAa = unname(varAas)))
  }
}

# [END]
