test_that("Handles invalid input", {
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  expect_error(processVariantData("/", protein = TRUE))
  expect_error(processVariantData(EvePath, protein = FALSE))
})

test_that("Results are of expected length", {
  varDataProtPath <- system.file("extdata", "variant_data_protein.csv", package = "variantMapper")
  processedVarDataProt <- processVariantData(varDataProtPath, protein = TRUE)
  expect_equal(nrow(processedVarDataProt), 7)
  expect_equal(ncol(processedVarDataProt), 3)

  varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
  processedVarDataGen <- processVariantData(varDataGenPath, protein = FALSE)
  expect_equal(nrow(processedVarDataGen), 7)
  expect_equal(ncol(processedVarDataGen), 5)
})

test_that("Results have expected columns", {
  protCols <- c("wtAa", "resPos", "varAa")
  genCols <- c("CHROM", "start", "end", "REF", "ALT")

  varDataProtPath <- system.file("extdata", "variant_data_protein.csv", package = "variantMapper")
  processedVarDataProt <- processVariantData(varDataProtPath, protein = TRUE)
  expect_equal(colnames(processedVarDataProt), protCols)

  varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
  processedVarDataGen <- processVariantData(varDataGenPath, protein = FALSE)
  expect_equal(colnames(processedVarDataGen), genCols)
})
