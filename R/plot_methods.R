plotupca <- function(obj, ...) {
  o <- out(obj)
  matplot(o[, 1], o[, -1], type = "l", ...)
  legend("topright", legend = c("u", "v", "w"), lty = 1:3,, bg = "white",
    col = 1:3)
}

#' @importFrom magrittr %>%
plotnp <- function(x, ...) UseMethod("plotnp")
plotnp.default <- function(x, ...) "Unknown class"
plotnp.odeModel <- function(obj, ...) {

  o <- tibble::as.tibble(out(obj)) %>%
    dplyr::select(time, N, P) %>%
    tidyr::gather(species, rho, N, P)

  cols <- c("N" = "darkgreen", "P" = "black")
  axis_labels <- labs(x = "Timestep", y = "Density",
    colour = "Species")

  g <- ggplot2::ggplot(o, aes(time, rho)) +
    geom_line(aes(colour = factor(species)))
  g + scale_colour_manual(values = cols) +
    axis_labels +
    ylim(0, 1)

}
plotnp.avg_scenarii <- function(data, x, threshold = 10^-3, debug_mode = FALSE, ...) {

  #var <- enquo(var)
  var <- quos(...)

  params <- data[["param"]]
  inits <- data[["inits"]]
  total_var  <- names(inits[[1]])

  run <- data[["run"]]

  #Replace low values by 0
  run %<>%
    tidyr::gather(species, rho, total_var) %>%
    dplyr::mutate(
      rho = replace(rho, rho < threshold, 0)
      ) %>%
    purrr::modify_at("species", as.factor)

  if (any(is.nan(run$rho))) {
    run %<>% dplyr::mutate(rho = replace(rho, is.nan(rho), 0))
    warning("NaNs produced in simulations have been replaced by 0")
  }

  run %<>% tidyr::spread(species, rho) %>%
    tidyr::gather(species, rho, !!! var) %>%
    purrr::modify_at("species", as.factor)

  g <- ggplot2::ggplot(run,
    aes_(y = ~rho, x = substitute(x), color = ~species, linetype = ~scenario)) +
  ggplot2::scale_linetype_manual(values = c("dotted", "solid")) +
  ggplot2::scale_colour_manual(values = c("green", "black")) +
  ggplot2::ylim(0,1) +
  theme_diagram() + # Provide a bifurcation theme 
  ggplot2::geom_line(alpha = .15, size = 3)

  if (debug_mode) {
    return(run)
  } else {
    return(g)
  }

}
plotnp.scenarii <- function(data, threshold = 10^-3, debug_mode = FALSE, alpha = .45) {
  #TODO: generalise for other parameters than b

  params <- data[["param"]]
  data %<>% .[["run"]]

  #Replace low values by 0
  var_to_drop <- c("NP", "NN", "PP", "status")
  data %<>%
    purrr::modify_at(var_to_drop, ~NULL) %>%
    gather(species, rho, N, P) %>%
    mutate(
      rho = replace(rho, rho < threshold, 0)
      ) %>%
    purrr::modify_at(c("scenario", "species"), as.factor)

  if (any(is.nan(data$rho))) {
    data %<>% mutate(rho = replace(rho, is.nan(rho), 0))
    warning("NaNs produced in simulations have been replaced by 0")
  }

  g <- ggplot(data,
    aes(y = rho, x = b, color = species, linetype = scenario)) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_colour_manual(values = c("green", "black")) +
  geom_line(alpha = alpha, size = 3) +
  theme_diagram() +
  ylim(0, 1) +
  facet_wrap(~ g)


  if (debug_mode) {
    return(data)
  } else {
    return(g)
  }

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
plot_diagram <- function(x, ...) UseMethod("plot_diagram")
plot_diagram.default <- function(x) "Unknown class"
plot_diagram.states <- function(
  data,
  param = c(x = "gamma1", y = "g"),
  possible_states = c("coexistence", "nurse", "protégée", "extinct", "warning"),
  col_states = c("orange", "green", "black", "yellow", "grey"),
  debug_mode = FALSE, ...) {

  cols <- col_states
  names(cols) <- possible_states

  g <- ggplot2::ggplot(data,
    aes_string(x = param["x"], y = param["y"], fill = "state")) +
  ggplot2::geom_raster() +
  theme_diagram() +
  ggplot2::scale_fill_manual(
    values = col_states,
    limits = possible_states
    )

  return(g)
}
plot_diagram.gradient <- function (
  data,
  param = c(x = "b", y = "g"),
  possible_states = c("coexistence", "nurse", "protegee", "desert", "warning"),
  col_states = c("orange", "green", "black", "yellow", "grey"),
  debug_mode = FALSE, ...) {

  params <- data[["param"]]
  data <- compute_states(data, param, type = "single")

  cols <- col_states
  names(cols) <- possible_states

  g <- ggplot2::ggplot(data$run,
    aes_string(x = param["x"], y = param["y"], fill = "state")) +
  ggplot2::geom_raster() +
  theme_diagram() +
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
plot_diagram.scenarii <- function (
  data,
  param = c(x = "b", y = "g"),
  possible_states = c("coexistence", "nurse", "protegee", "desert", "warning"),
  col_states = c(coexistence = "orange", nurse = "green", protegee = "black",
  desert = "#C19A6B", unkown = "white"),
  debug_mode = FALSE, fill = "single", ...) {
  if (fill %in% c("single", "double_states")) {
    data <- compute_states(data,
      param = param,
      fill = fill)

    g <- ggplot2::ggplot(data$run,
      aes_string(x = param["x"], y = param["y"], fill = "state")) +
    ggplot2::geom_raster() +
    theme_diagram() +
    ggplot2::scale_fill_manual(
      values = col_states,
      limits = possible_states
      )

    if (fill == "single"){
     g <- g + facet_grid(. ~ scenario)
    }
  } else {
    if (fill == "cnp") {
      data %<>% compute_occurences(.)
    }

    g <- ggplot2::ggplot(data$run,
      aes_string(x = param["x"], y = param["y"], fill = fill)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient2(
      low = "blue", high = "red", limits = c(0, 1.5), midpoint = 1
      ) +
    theme_diagram()
  }
  if (debug_mode) {
    return(data)
  }
  return(g)
}
theme_diagram <- function(base_size = 12, base_family = "Helvetica"){
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      #line = element_line(colour="black"),
      #text = element_text(colour="black"),
      axis.title = element_text(size = base_size),
      #axis.text = element_text(colour="black", size=8),
      #strip.text = element_text(size=12),
      legend.key=element_rect(colour=NA, fill =NA),
      panel.grid = element_blank(),   
      panel.border = element_rect(fill = NA, colour = "black", size = 1),
      panel.background = element_rect(fill = "white", colour = "black"),
      strip.background = element_rect(fill = NA)
      )
}
