# helper function to get unique wtAa and the resPos
uniqueWtAaPos <- function(eveData) {
  # Note- wtAa and resPos are columns in eveData
  # TODO: Include check to make sure they are in the eveData
  wtAaPos <- dplyr::distinct(eveData, wtAa, resPos)
  return(wtAaPos)
}

# helper function to find out which nt is mutated (1, 2, or 3)
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

# helper function to construct alternate codon sequence
constructAltSeq <- function(refData, ntPos, genomicCoord, varNt) {
  # TODO: POS is column in refData, check if it is present
  refInfo <- dplyr::filter(refData, POS == genomicCoord) # get the matching row
  refCodon <- refInfo[1,]$REF
  resPos <- refInfo[1,]$resPos

  # produce the alternative codon seq
  altCodon <- refInfo[1,]$REF # initially assign it ref codon seq
  stringr::str_sub(altCodon, ntPos, ntPos) <- varNt # change to alternate codon
                                                    # seq
  return(tibble(refCodon = refCodon, altCodon = altCodon, resPos = resPos))
}

scoreVariants <- function(eveData, variantData, protein = TRUE) {

  # eveData is the tibble
  # return a vector with the EVE score for all the proteins with EVE score

  # call helper to get pairs of distinct wtAa and residue positions
  wtAaPos <- uniqueWtAaPos(eveData)
  numResidues <- nrow(wtAaPos)
  eveScores <- vector("numeric", numResidues) # vector to return

  # those with no mutation i.e. wildtype have EVE score 0 (benign), and those
  # variants who aren't scored are assigned NaN

  if (isTRUE(protein)) {
    # no translation required

    # TODO: do a null check for variants
    variants <- dplyr::inner_join(eveData, variantData) # variants with EVE data

    eveCol <- which(colnames(variants) == "EVE") # get column with EVE score

    # loop through the variant data
    for (i in 1:nrow(variantData)) {
      mut <- variantData[i,]

      # get a subset of all the rows in EVE data with that wtAa, resPos, varAa
      varSubset <- dplyr::filter(variants,
                                 wtAa == mut$wtAa,
                                 resPos == mut$resPos,
                                 varAa == mut$varAa)

      # score the mutation
      if (nrow(varSubset) == 0) {
        # variant isn't scored
        eveScores[mut$resPos] <- NaN
      } else if (nrow(varSubset) == 1) {
        eveScores[mut$resPos] <- varSubset[1, eveCol]
      } else {
        # assign EVE score as the average
        # TODO: look into if all variants mutate to same aa but different codon
        # have same EVE score
        eveScores[mut$resPos] <- mean(varSubset$EVE)
      }
    }
  } else {
    # if it is in genomic coordinates
  }
  return(eveScores)

}
