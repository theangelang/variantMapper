% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoreVariant.R
\name{scoreVariant}
\alias{scoreVariant}
\title{Scores a variant of a protein based on weighted average}
\usage{
scoreVariant(
  eveScores,
  posWeights = rep(1/length(eveScores), length(eveScores))
)
}
\arguments{
\item{eveScores}{A vector of doubles containing EVE scores for each residue
position that has a possible EVE score.}

\item{posWeights}{A vector of doubles representing the weights to give each
residue position when calculating the weighted average of EVE scores for the
protein variant.  The default is using equal weights for each residue
position i.e. 1/(number of residues with EVE scores).  This is used when an
argument is not passed into the posWeights parameter.}
}
\value{
A double as the weighted EVE score for the protein variant.
}
\description{
A function that returns the weighted average EVE score of a protein variant.
}
\details{
EVE (Evolutionary model of Variant Effect) is an unsupervised machine
learning model shown to be accurate in predicting pathogenicity of missense
variants.  It uses multiple sequence alignments and doesn't rely on knowledge
of protein function to do so.  More about EVE can be found here
(https://evemodel.org/).

The EVE score assigned is continuous on the interval zero to one.  An EVE
score of zero indicates benign while an EVE score of one is most pathogenic.
}
\examples{
# Examples:
# First process the EVE data and variant data.
EvePath <- system.file("extdata", "NRX1B_HUMAN_SUBSET.vcf", package = "variantMapper")
EveData <- processEveData(EvePath)
EveData

# If the data is in protein form.
varDataProtPath <- system.file("extdata", "variant_data_protein.csv", package = "variantMapper")
varDataProt <- processVariantData(varDataProtPath, protein = TRUE)
varDataProt

eveScoresProt <- getEveScores(EveData, varDataProt, protein = TRUE)
eveScoresProt

# Get the weighted average EVE score to see overall how pathogenic this variant is
avgEveScore <- scoreVariant(eveScoresProt$eveScores) # using the default of uniform weights
avgEveScore

# If the data is in genomic form.
varDataGenPath <- system.file("extdata", "variant_data_genomic.csv", package = "variantMapper")
varDataGen <- processVariantData(varDataGenPath, protein = FALSE)
varDataGen

eveScoresGen <- getEveScores(EveData, varDataGen, protein = FALSE)
eveScoresGen

# Get the weighted average EVE score to see overall how pathogenic this variant is
# Only assign weight to the first few residue positions
weightsUnequal <- rep(0, length(eveScoresGen$eveScores))
weightsUnequal[1] <- 0.8
weightsUnequal[2] <- 0.1
weightsUnequal[3] <- 0.1
weightedAvgEveScore <- scoreVariant(eveScoresGen$eveScores, weightsUnequal)
weightedAvgEveScore

}
\references{
1. R Core Team (2022). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.
URL https://www.R-project.org/.

2. Frazer, J. et al. (2021). Disease variant prediction with deep generative
models of evolutionary data. \emph{Nature. 599}. 91-95.
}
