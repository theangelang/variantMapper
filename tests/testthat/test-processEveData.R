test_that("Handles invalid input", {
  variantPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
  expect_error(processEveData("/"))
  expect_error(processEveData(variantPath))
})

test_that("Results are of expected length", {
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  processedEveData <- processEveData(EvePath)
  expect_equal(1824, nrow(processedEveData))
  expect_equal(26, ncol(processedEveData))
})

test_that("Results have expected columns", {
  columns <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "Key", "EVE",
               "EnsTranscript", "RevStr", "ProtMut", "Class10", "Class20", "Class25",
               "Class30", "Class40", "Class50", "Class60", "Class70", "Class75",
               "Class80", "Class90", "wtAa", "resPos", "varAa")
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  processedEveData <- processEveData(EvePath)
  expect_equal(columns, colnames(processedEveData))
})
