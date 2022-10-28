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
#' @import vcfR dplyr tibble
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

  return(eveData)

}

#[END]