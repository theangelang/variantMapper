#' Launch Shiny App for variantMapper
#'
#' A function that launches the Shiny app for variantMapper, the code for the
#' Shiny app haas been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value but will open up a Shiny page.
#'
#' @examples
#' \dontrun{
#' variantMapper::runVariantMapper()
#' }
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
