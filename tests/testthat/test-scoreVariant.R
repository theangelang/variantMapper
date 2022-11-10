test_that("Handles invalid input", {
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  eveData <- processEveData(EvePath)

  varDataProtPath <- system.file("extdata", "variant_data_protein.csv", package = "variantMapper")
  varDataProt <- processVariantData(varDataProtPath, protein = TRUE)

  varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
  varDataGen <- processVariantData(varDataGenPath, protein = FALSE)

  eveScoresProt <- getEveScores(eveData, varDataProt, protein = TRUE)
  eveScoresGen <- getEveScores(eveData, varDataGen, protein = FALSE)

  shorterPosWeights <- c(1/(length(eveScoresProt) - 1), length(eveScoresProt) - 1)
  expect_warning(scoreVariant(eveScoresProt$eveScores, shorterPosWeights))

  invalidPosWeights <- c(1, length(eveScoresGen))
  expect_warning(scoreVariant(eveScoresGen$eveScores, invalidPosWeights))

  invalidEveScores <- c(NaN, length(eveScoresProt))
  expect_warning(scoreVariant(eveScoresGen$eveScores, invalidEveScores))
})

test_that("Results are of expected type", {
  EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
  eveData <- processEveData(EvePath)

  varDataProtPath <- system.file("extdata", "variant_data_protein.csv", package = "variantMapper")
  varDataProt <- processVariantData(varDataProtPath, protein = TRUE)

  varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
  varDataGen <- processVariantData(varDataGenPath, protein = FALSE)

  eveScoresProt <- getEveScores(eveData, varDataProt, protein = TRUE)
  eveScoresGen <- getEveScores(eveData, varDataGen, protein = FALSE)

  weightsUnequal <- c(0.3, 0.2, 0.5)

  expect_type(scoreVariant(eveScoresProt$eveScores), "double")
  expect_type(scoreVariant(eveScoresGen$eveScores[1:3], weightsUnequal), "double")
})
