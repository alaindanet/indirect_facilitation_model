
plotupca <- function(obj, ...) {
  o <- out(obj)
  matplot(o[, 1], o[, -1], type = "l", ...)
  legend("topright", legend = c("u", "v", "w"), lty = 1:3,, bg = "white",
    col = 1:3)
}


#' @importFrom magrittr %>%
plotnp <- function(obj, ...) {
  o <- tibble::as.tibble(out(obj)) %>%
    dplyr::select(time, N, P) %>%
    tidyr::gather(species, rho, N, P)

  cols <- c("N" = "darkgreen", "P" = "black")
  axis_labels <- labs(x = "Timestep", y = "Density",
    colour = "Species")

  g <- ggplot2::ggplot(o, aes(time, rho)) +
    geom_line(aes(colour = factor(species)))
  g + scale_colour_manual(values = cols) +
    axis_labels

}
