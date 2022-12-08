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

test_that("Helper function checkInputData results are of expected type", {
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  processedEveData <- processEveData(EvePath)

  varDataProtPath <- system.file("extdata", "variant_data_protein.csv", package = "variantMapper")
  varDataProt <- processVariantData(varDataProtPath, protein = TRUE)

  varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
  varDataGen <- processVariantData(varDataGenPath, protein = FALSE)

  trueResultEve <- checkInputData(colnames(processedEveData), "EVE")
  trueResultVarProt <- checkInputData(colnames(varDataProt), "protein")
  trueResultVarGen <- checkInputData(colnames(varDataGen), "genomic")

  falseResultEve <- checkInputData(colnames(processedEveData), "protein")
  falseResultVarProt <- checkInputData(colnames(varDataProt), "genomic")
  falseResultVarGen <- checkInputData(colnames(varDataGen), "EVE")

  expect_true(trueResultEve)
  expect_true(trueResultVarProt)
  expect_true(trueResultVarGen)

  expect_false(falseResultEve)
  expect_false(falseResultVarProt)
  expect_false(falseResultVarGen)

  expect_error(checkInputData(colnames(processedEveData), "other"))
})

test_that("Helper function uniqueWtAa results are of expected length and columns", {
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  processedEveData <- processEveData(EvePath)

  columns <- c("wtAa", "resPos")

  result <- uniqueWtAaPos(processedEveData)

  expect_equal(nrow(result), 32)
  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), columns)
})

test_that("Helper function findVariantPosition results are of expected length", {
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  processedEveData <- processEveData(EvePath)

  varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
  varDataGen <- processVariantData(varDataGenPath, protein = FALSE)

  coordinates <- dplyr::distinct(processedEveData, POS)

  result <- findVariantPosition(coordinates, varDataGen[[1,"start"]])

  expect_equal(length(result), 2)
})

test_that("Helper function constructAltSeq results are of expected length and value", {
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  processedEveData <- processEveData(EvePath)

  varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
  varDataGen <- processVariantData(varDataGenPath, protein = FALSE)

  columns <- c("refCodon", "altCodon", "resPos")

  refData <- dplyr::distinct(processedEveData, POS, REF, resPos)

  result1 <- constructAltSeq(refData = refData,
                             ntPos = 2,
                             genomicCoord = 50346854,
                             varNt = "T")
  result2 <- constructAltSeq(refData = refData,
                             ntPos = 3,
                             genomicCoord = 50346851,
                             varNt = "T")
  result3 <- constructAltSeq(refData = refData,
                             ntPos = 1,
                             genomicCoord = 50346845,
                             varNt = "G")

  expect_equal(ncol(result1), 3)
  expect_equal(nrow(result1), 1)
  expect_equal(colnames(result1), columns)
  expect_equal(result1[[1, "altCodon"]], "GTC")

  expect_equal(ncol(result2), 3)
  expect_equal(nrow(result2), 1)
  expect_equal(colnames(result2), columns)
  expect_equal(result2[[1, "altCodon"]], "CAT")

  expect_equal(ncol(result3), 3)
  expect_equal(nrow(result3), 1)
  expect_equal(colnames(result3), columns)
  expect_equal(result3[[1, "altCodon"]], "GCA")
})
