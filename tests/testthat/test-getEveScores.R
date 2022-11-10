test_that("Handles invalid input", {
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  eveData <- processEveData(EvePath)

  varDataProtPath <- system.file("extdata", "variant_data_protein.csv", package = "variantMapper")
  varDataProt <- processVariantData(varDataProtPath, protein = TRUE)

  varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
  varDataGen <- processVariantData(varDataGenPath, protein = FALSE)

  expect_error(getEveScores(eveData[,1:3], varDataProt, protein = TRUE))
  expect_error(getEveScores(eveData, varDataProt[,1:2], protein = TRUE))
  expect_error(getEveScores(eveData, varDataGen[,1:3], protein = FALSE))
  expect_error(getEveScores(eveData, varDataProt, protein = FALSE))
})

test_that("Results are of expected length", {
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  eveData <- processEveData(EvePath)

  varDataProtPath <- system.file("extdata", "variant_data_protein.csv", package = "variantMapper")
  varDataProt <- processVariantData(varDataProtPath, protein = TRUE)

  varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
  varDataGen <- processVariantData(varDataGenPath, protein = FALSE)

  expect_equal(nrow(getEveScores(eveData, varDataProt, protein = TRUE)), 32)
  expect_equal(ncol(getEveScores(eveData, varDataProt, protein = TRUE)), 4)

  expect_equal(nrow(getEveScores(eveData, varDataGen, protein = FALSE)), 32)
  expect_equal(ncol(getEveScores(eveData, varDataGen, protein = FALSE)), 4)
})

test_that("Results have expected columns", {
  columns <- c("eveScores", "resPos", "wtAa", "varAa")

  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  eveData <- processEveData(EvePath)

  varDataProtPath <- system.file("extdata", "variant_data_protein.csv", package = "variantMapper")
  varDataProt <- processVariantData(varDataProtPath, protein = TRUE)

  varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
  varDataGen <- processVariantData(varDataGenPath, protein = FALSE)

  expect_equal(colnames(getEveScores(eveData, varDataProt, protein = TRUE)), columns)
  expect_equal(colnames(getEveScores(eveData, varDataGen, protein = FALSE)), columns)
})
