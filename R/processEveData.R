#' Processes the data from EVE in vcf format
#'
#' A function that checks if the data provided is in the right format (vcf) and
#' processes it for future use.
#'
#' @param filePath The filepath to the vcf file from EVE containing information
#' about the gene of interest such as the EVE scores.  This can be obtained
#' from the EVE website by searching for the gene then selecting vcf under the
#' download section.
#'
#' @return Returns a tibble with the data from the vcf file provided.
#'
#' @import vcfR dplyr tibble stringr
processEveData <- function(filePath) {

  # TODO: convert the column types to integers where appropriate
  if (!file.exists(filePath)) {
    stop("File not found.  Please provide a valid filepath for the EVE data.")
  }

  if (!grepl("*\\.vcf$",filePath)){
    stop("File format not valid, please provide a .vcf file.")
  }

  # read in the vcf file
  readData <- vcfR::read.vcfR(filePath)

  # get information for chromosome, genomic position, reference seq, mutated seq
  seqInfo <- vcfR::getFIX(readData)
  seqInfoTibble <- tibble::as_tibble(seqInfo) # convert to tibble

  # extract information from INFO column of vcf file
  eveScores <- vcfR::extract_info_tidy(readData)

  # combine 2 tibbles to 1
  eveData <- dplyr::bind_cols(seqInfoTibble, eveScores)

  # convert the CHROM and POS columns to numeric
  i <- c(1, 2)
  eveData[ , i] <- apply(eveData[ , i],
                         2,
                         function(x) as.numeric(as.character((x))))

  # convert the RevStr to type logical
  eveData$RevStr <- as.logical(as.character(eveData$RevStr))

  protChange <- getProtMutInfo(eveData$ProtMut)

  eveData <- dplyr::bind_cols(eveData, protChange)

  return(eveData)

}

# helper function to extract the wild type, residue position, and variant amino
# acid from the ProtMut column
getProtMutInfo <- function(protMut) {

  aaChange <- sapply(strsplit(protMut, "_"), "[", 2) # keep everything after _

  # get the wild type amino acid
  wtAa <- substr(aaChange, 1, 1)

  # get the variant amino acid
  varAa <- stringr::str_sub(aaChange, -1) # from stringr package

  resPos <- gsub("[^0-9.-]", "", aaChange) # replace all the characters with ""
  resPos <- as.numeric(as.character(resPos))# convert to numeric

  res <- dplyr::tibble(wtAa, resPos, varAa)

  return(res)
}

#[END]
