#' Processes the data from EVE in vcf format
#'
#' A function that checks if the data provided is in the right format (vcf) and
#' processes it for future use.
#'
#' This data (vcf file) should come from the EVE website
#' (https://evemodel.org/download/protein) under the 'Download' tab.
#' On this page enter in the gene of interest and select VCFs.
#'
#' EVE (Evolutionary model of Variant Effect) is an unsupervised machine
#' learning model shown to be accurate in predicting pathogenicity of missense
#' variants.  It uses multiple sequence alignments and doesn't rely on knowledge
#' of protein function to do so.  More about EVE can be found here
#' (https://evemodel.org/).
#'
#' The EVE score assigned is continuous on the interval zero to one.  An EVE
#' score of zero indicates benign while an EVE score of one is most pathogenic.
#'
#' @param filePath A character vector of the filepath to the vcf file from EVE
#' containing information about the gene of interest such as the EVE scores.
#' This can be obtained from the EVE website
#' (https://evemodel.org/download/protein) by searching for the gene then
#' selecting vcf under the download section.
#'
#' @return A tibble with the EVE data from the vcf file provided.
#'
#' @export
#'
#' @examples
#' # Example:
#' # Get the EVE data in the correct format by accessing raw EVE data provided
#' # for a gene
#'
#' EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
#' EveData <- processEveData(EvePath)
#' EveData
#'
#' @references
#' 1. R Core Team (2022). R: A language and environment for statistical computing.
#' R Foundation for Statistical Computing, Vienna, Austria.
#' URL https://www.R-project.org/.
#'
#' 2. Knaus BJ, Grünwald NJ (2017). “VCFR: a package to manipulate and visualize
#' variant call format data in R.” \emph{Molecular Ecology Resources, 17}(1),
#' 44-53. ISSN 757, http://dx.doi.org/10.1111/1755-0998.12549.
#'
#' 3. Müller K, Wickham H (2022). \emph{tibble: Simple Data Frames}. R package
#' version 3.1.8, https://CRAN.R-project.org/package=tibble.
#'
#' 4. Wickham H (2022). \emph{stringr: Simple, Consistent Wrappers for Common
#' String Operations}. R package version 1.4.1,
#' <https://CRAN.R-project.org/package=stringr>.
#'
#' 5. Frazer, J. et al. (2021). Disease variant prediction with deep generative
#' models of evolutionary data. \emph{Nature. 599}. 91-95.
#'
#' @import vcfR dplyr tibble stringr

processEveData <- function(filePath) {

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

  # convert the POS columns to numeric
  eveData[ , 2] <- apply(eveData[ , 2],
                         2,
                         function(x) as.numeric(as.character((x))))

  # convert the RevStr to type logical
  eveData$RevStr <- as.logical(as.character(eveData$RevStr))

  protChange <- getProtMutInfo(eveData$ProtMut)

  eveData <- dplyr::bind_cols(eveData, protChange)

  return(eveData)

}

#' Helper function to extract wild type, residue position, and variant amino
#' acid from ProtMut column in the EVE dataset.
#'
#' @param protMut A character vector of form
#' UniprotId_(WildTypeAminoAcid)(ResiduePosition)(VariantAminoAcid) without the
#' brackets.
#'
#' @return A tibble containing the wild type amino acid, residue
#' position, and variant amino acid.
#'
#' @noRd

getProtMutInfo <- function(protMut) {

  aaChange <- sapply(strsplit(protMut, "_"), "[", 2) # keep everything after _

  # get the wild type amino acid
  wtAa <- substr(aaChange, 1, 1)

  # get the variant amino acid
  varAa <- stringr::str_sub(aaChange, -1)

  resPos <- gsub("[^0-9.-]", "", aaChange) # replace all the characters with ""
  resPos <- as.numeric(as.character(resPos)) # convert to numeric

  res <- dplyr::tibble(wtAa, resPos, varAa)

  return(res)
}

#[END]
