#' Processes variant data in csv format
#'
#' A function that checks if the data provided is in the right format (csv), if
#' it has the columns required for downstream analysis depending on if
#' the variant data is in genomic or protein coordinates, and processes it for
#' future use.
#'
#' @param filePath The filepath to the csv file containing variant information
#' about the variant and gene of interest.
#'
#' @param protein This argument specifies whether the input data is in genomic
#' or protein coordinates.
#'
#' @return Returns a dataframe with the data from the csv file provided.
#'
#' @importFrom utils read.csv

processVariantData <- function(filePath, protein = TRUE) {

  if (!file.exists(filePath)) {
    stop("File not found.  Please provide a valid filepath for the variant
         data.")
  }

  if (!grepl("*\\.csv$",filePath)){
    stop("File format not valid, please provide a .csv file.")
  }

  # read in variant data
  varData <- read.csv(file = filePath)

  cols <- colnames(varData)

  # check if data is in protein coordinates
  if (isTRUE(protein)) {

    if ("wt_aa" %in% cols & "position" %in% cols & "mt_aa" %in% cols) {
      return (varData)
    } else {
      stop("File is missing columns, please provide a .csv file with columns
           'wt_aa', 'position', and 'mt_aa.")
    }
  } else {

    if ("chrom" %in% cols & "start" %in% cols & "end" %in% cols &
        "ref_allele" %in% cols & "alt_allele" %in% cols & "var_type" %in% cols) {
      return (varData)
    } else {
      stop("File is missing columns, please provide a .csv file with columns
           'chrom', 'start', 'end', 'ref_allele', 'alt_allele', and
           'var_type'.")
    }
  }
}

# [END]
