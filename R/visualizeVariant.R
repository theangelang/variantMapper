#' Visualizes EVE score for a protein variant.
#'
#' @import ggplot2

visualizeVariant <- function(eveInfo) {
  # TODO: add check to make sure has those columns
  p <- ggplot(eveInfo, aes(x=resPos, y=eveScores)) +
    geom_segment( aes(x=resPos, xend=resPos, y=0, yend=eveScores), color="grey") +
    geom_point(aes(color=eveScores), size=4) +
    scale_colour_gradient2(
      low = "steelblue1",
      mid = "gray",
      high = "firebrick1",
      midpoint = 0.5) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(x = "Residue Position", y = "EVE Score", color = "EVE Score")
  return(p)
}
