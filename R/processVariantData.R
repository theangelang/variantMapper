#' Processes variant data in csv format
#'
#' A function that checks if the data provided is in the right format (csv), if
#' it has the columns required for downstream analysis depending on if
#' the variant data is in genomic or protein coordinates, and processes it for
#' future use.
#'
#' The data (csv file) should represent single nucleotide variants.
#'
#' If the single nucleotide variant (SNV) data is provided in protein form it
#' should have three columns.The first being 'wt_aa' with the single amino acid
#' code for the wild type amino acid.  The second is 'position' containing
#' integers representing residue position.  The third column is 'mt_aa'which is
#' the amino acid variant at this particular residue position.  It also contains
#' the one letter amino acid code.  For an example please see
#' 'variant_data_protein.csv' at
#' https://github.com/theangelang/variantMapper/tree/master/inst/extdata
#'
#' If the single nucleotide variant data is in genomic form it should be a csv
#' file with five columns.  The first being 'chrom' which are either integers or
#' X or Y representing chromosome location.  The second being 'start' which has
#' integers indicating start position of the SNV.  The third is 'end',
#' containing integers representing the end position of the SNV.  The fourth
#' column is 'ref_allele' which has characters representing the wildtype
#' nucleotide.  The fifth column is 'alt_allele' which has characters
#' representing the variant nucleotide.  For an example please see
#' 'variant_data_genomic.csv' at
#' https://github.com/theangelang/variantMapper/tree/master/inst/extdata"
#'
#' This package aims to make variant pathogenicity predictions more accessible
#' and it does that through annotating variants with an EVE score. In order to
#' annotate a variant of interest with EVE scores this function must
#' be used prior to the 'getEveScores' function to properly process the variant
#' data.
#'
#' EVE (Evolutionary model of Variant Effect) is an unsupervised machine
#' learning model shown to be accurate in predicting pathogenicity of missense
#' variants.  It uses multiple sequence alignments and doesn't rely on knowledge
#' of protein function to do so.  More about EVE can be found here
#' (https://evemodel.org/).
#'
#' @param filePath A character vector of the filepath to the csv file containing
#' variant information for the gene of interest.
#'
#' @param protein A logical value specifying whether the input data is in
#' genomic or protein coordinates.  If TRUE the data contains protein
#' information and protein coordinates; this is the default.  If FALSE it
#' contains genomic information and genomic coordinates.
#'
#' @return A tibble with the data from the csv file provided.
#'
#' @export
#'
#' @examples
#' # Example
#' # Get the variant data in the right format depending on if it is protein or
#' # genomic data.
#'
#' # If the data is in protein form.
#' varDataProtPath <- system.file("extdata", "variant_data_protein.csv", package = "variantMapper")
#' varDataProt <- processVariantData(varDataProtPath, protein = TRUE)
#' varDataProt
#'
#' # If the data is in genomic form.
#' varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
#' varDataGen <- processVariantData(varDataGenPath, protein = FALSE)
#' varDataGen
#'
#' @references
#' 1. R Core Team (2022). R: A language and environment for statistical
#' computing. R Foundation for Statistical Computing, Vienna, Austria.
#' URL https://www.R-project.org/.
#'
#' 2. Müller K, Wickham H (2022). \emph{tibble: Simple Data Frames}. R package
#' version 3.1.8, https://CRAN.R-project.org/package=tibble.
#'
#' 3. Frazer, J. et al. (2021). Disease variant prediction with deep generative
#' models of evolutionary data. \emph{Nature. 599}. 91-95.
#'
#' @import tibble
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
  varData <- utils::read.csv(file = filePath)

  cols <- colnames(varData)

  # check if data is in protein coordinates
  if (isTRUE(protein)) {
    if ("wt_aa" %in% cols & "position" %in% cols & "mt_aa" %in% cols) {
      varDataTibble <- tibble::tibble(wtAa = varData$wt_aa,
                                      resPos = varData$position,
                                      varAa = varData$mt_aa)
      return (varDataTibble)

    } else {
      stop("File is missing columns, please provide a .csv file with columns
           'wt_aa', 'position', and 'mt_aa.")
    }
  } else {
    if ("chrom" %in% cols & "start" %in% cols & "end" %in% cols &
        "ref_allele" %in% cols & "alt_allele" %in% cols) {
      varDataTibble <- tibble::tibble(CHROM = varData$chrom,
                              start = varData$start,
                              end = varData$end,
                              REF = varData$ref_allele,
                              ALT = varData$alt_allele)

      if (!identical(varDataTibble[["start"]], varDataTibble[["end"]])) {
        warning("The variant data contains variants other than single nucleotide
                variants which can lead to unusual results.")
      }

      return (varDataTibble)

    } else {
      stop("File is missing columns, please provide a .csv file with columns
           'chrom', 'start', 'end', 'ref_allele', 'alt_allele', and
           'var_type'.")
    }
  }
}

# [END]
