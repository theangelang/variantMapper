test_that("Handles invalid input", {
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  eveData <- processEveData(EvePath)

  varDataProtPath <- system.file("extdata", "variant_data_protein.csv", package = "variantMapper")
  varDataProt <- processVariantData(varDataProtPath, protein = TRUE)

  varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
  varDataGen <- processVariantData(varDataGenPath, protein = FALSE)

  eveScoresProt <- getEveScores(eveData, varDataProt, protein = TRUE)
  eveScoresGen <- getEveScores(eveData, varDataGen, protein = FALSE)

  expect_error(visualizeVariant2(eveScoresProt[,1:3], eveScoresGen[,1:3]))
  expect_error(visualizeVariant2(eveScoresProt, eveScoresGen[,1:3], geneName = "NRXN1", aboveZeroOnly = TRUE))
  expect_error(visualizeVariant2(eveScoresProt[,1:3], eveScoresGen, geneName = "NRXN1"))
})
