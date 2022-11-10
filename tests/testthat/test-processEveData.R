test_that("Handles invalid input", {
  variantPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
  expect_error(processEveData("/"))
  expect_error(processEveData(variantPath))
})

test_that("Results are of expected length", {
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  processedEveData <- processEveData(EvePath)
  expect_equal(nrow(processedEveData), 1824)
  expect_equal(ncol(processedEveData), 26)
})

test_that("Results have expected columns", {
  columns <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "Key", "EVE",
               "EnsTranscript", "RevStr", "ProtMut", "Class10", "Class20", "Class25",
               "Class30", "Class40", "Class50", "Class60", "Class70", "Class75",
               "Class80", "Class90", "wtAa", "resPos", "varAa")
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  processedEveData <- processEveData(EvePath)
  expect_equal(colnames(processedEveData), columns)
})

test_that("Helper function results are of expected length and columns", {
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  processedEveData <- processEveData(EvePath)

  columns <- c("wtAa", "resPos", "varAa")

  result <- getProtMutInfo(processedEveData$ProtMut)

  expect_equal(nrow(result), 1824)
  expect_equal(ncol(result), 3)
  expect_equal(colnames(result), columns)
})

