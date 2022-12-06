#' Launch Shiny App for variantMapper
#'
#' A function that launches the Shiny app for variantMapper, the code for the
#' Shiny app has been placed in \code{./inst/shiny-scripts}.
#'
#' The Shiny app launched by this function allows for the use of this package's
#' functionality in a web interface.  In this Shiny app you will be able to
#' upload your variant data, assign variants of interest an EVE score, visualize
#' and compare variants.
#'
#' EVE (Evolutionary model of Variant Effect) is an unsupervised machine
#' learning model shown to be accurate in predicting pathogenicity of missense
#' variants.  It uses multiple sequence alignments and doesn't rely on knowledge
#' of protein function to do so.  More about EVE can be found here
#' (https://evemodel.org/).
#'
#' @return No return value but will open up a Shiny page.
#'
#' @examples
#' \dontrun{
#' variantMapper::runVariantMapper()
#' }
#'
#' @references
#' 1. R Core Team (2022). R: A language and environment for statistical
#' computing. R Foundation for Statistical Computing, Vienna, Austria.
#' URL https://www.R-project.org/.
#'
#' 2. Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Xie Y, Allen J,
#' McPherson J, Dipert A, Borges B (2022). \emph{shiny: Web Application
#' Framework for R}. R package version 1.7.3,
#' https://CRAN.R-project.org/package=shiny.
#'
#' 3. \emph{File Upload}. (2022). RStudio Shiny Gallery,
#' https://shiny.rstudio.com/gallery/file-upload.html.
#'
#' 4. \emph{Tabsets}. (2022). RStudio Shiny Gallery,
#' https://shiny.rstudio.com/gallery/tabsets.html.
#'
#' 5. \emph{Function reference}. (2022). RStudio Shiny,
#' https://shiny.rstudio.com/reference/shiny/1.7.3/.
#'
#' 6. Frazer, J. et al. Disease variant prediction with deep generative models
#' of evolutionary data. Nature. 599. 91-95 (2021).
#'
#' @export
#' @importFrom shiny runApp
#'

runVariantMapper <- function() {
  appDir <- system.file("shiny-scripts", package = "variantMapper")
  runShinyApp <- shiny::runApp(appDir, display.mode = "normal")
  return(runShinyApp)
}

# [END]
