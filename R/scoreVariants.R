# helper function to get unique wtAa and the resPos
uniqueWtAaPos <- function(eveData) {
  # Note- wtAa and resPos are columns in eveData
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
