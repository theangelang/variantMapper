#' Visualizes EVE score for a protein variant.
#'
#' @import ggplot2

visualizeVariant <- function(eveInfo) {
  # TODO: add check to make sure has those columns
  # TODO: handle the NaN values and turn to 0
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

#' Function to visualize EVE scores of two variants at once
#'
#' @import ggplot2

visualizeVariant2 <- function(eveInfo1, eveInfo2) {
  # this one compares variants from different samples onto one gene
  p <- ggplot(eveInfo1, aes(x=resPos, y=eveScores)) +
    geom_segment( aes(x=resPos, xend=resPos, y=0, yend=eveScores), color="grey") +
    geom_point(data=eveInfo1, aes(color="Variant 1"), size=4) +
    geom_point(data=eveInfo2, aes(color="Variant 2"), size=4) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(x = "Residue Position", y = "EVE Score", color = "EVE Score")
  return(p)
}
