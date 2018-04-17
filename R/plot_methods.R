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

#' Plot density over a gradient 
#'
#' @param data 
#' 
#' @return a plot 
#' @export
plotnp_gradient <- function(data, state_var = c("N", "P"), param = c("gamma1", "g"), ...) {

  var_to_drop <- names(data)[!(names(data) %in% c(state_var, param))]

  data %<>% purrr::modify_at(var_to_drop, ~NULL) %>%
    tidyr::gather(species, rho, state_var) %>%
    tidyr::gather(gradient, value, param)

  cols <- c("N" = "darkgreen", "P" = "black")

  g <- ggplot2::ggplot(data, aes(y = rho, x = value)) +
    geom_line(aes(colour = factor(species), linetype = factor(gradient))) +
    scale_colour_manual(values = cols)
  g


}

#' Plot diagram over a gradient 
#'
#' @param data 
#' @param possible_states chr vector specifying the different states   
#' @param col_states chr vector of the same length of the possible_states vector.
#' The order of the vector match the order of definition of the possible_states.
#' 
#' @return a plot
#' @export
plot_diagram <- function (
  data,
  param = c(x = "gamma1", y = "g"),
  possible_states = c("coexistence", "nurse", "protégée", "extinct", "warning"),
  col_states = c("orange", "green", "black", "yellow", "grey"),
  debug_mode = FALSE, ...) {

  var_to_drop <- names(data)[!(names(data) %in% c(param))]
  data %<>%
    dplyr::mutate(
      state = purrr::pmap_chr(
	list(nurse = N, protegee = P, sim_status = status),
	def_state)) %>%
    purrr::modify_at(var_to_drop, ~NULL) %>%
    dplyr::mutate(state = factor(state, levels = possible_states))

  cols <- col_states
  names(cols) <- possible_states

  g <- ggplot2::ggplot(data,
    aes_string(x = param["x"], y = param["y"], fill = "state")) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_manual(
    values = col_states,
    limits = possible_states
    )

  if (debug_mode) {
    return(data)
  } else {
    return(g)
  }

}
