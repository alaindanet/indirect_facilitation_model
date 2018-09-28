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
plot_diagram.default <- function(x, ...) "Unknown class"
plot_diagram.states_scenarii <- function (
  data,
  param = c(x = "b", y = "g"), debug_mode = FALSE) {

    g <- ggplot2::ggplot(data$run,
      aes_string(x = param["x"], y = param["y"], fill = "state")) +
    ggplot2::geom_raster() +
    theme_diagram() +
    ggplot2::scale_fill_manual(
      values = color_states()
      )

  if (debug_mode) {
    return(data)
  }
  return(g)
}
plot_diagram.avg_scenarii <- function (
  data,
  param = c(x = "b", y = "g"), debug_mode = FALSE, fill = "N") {

    g <- ggplot2::ggplot(data$run,
      aes_string(x = param["x"], y = param["y"], fill = fill)) +
    ggplot2::geom_raster() +
    theme_diagram()
  #+
    #ggplot2::scale_fill_manual(
      #values = color_states()
      #)

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

#' Color blind palette 
#'
#' Paul tol palette light qualitative, it can be found at 
#' https://personal.sron.nl/~pault/
#' 
#' @return a named vector of 9 colors
#' @export
compute_light_paul_palette <- function() {

  c(
    light_blue = "#77AADD",
    light_cyan = "#99DDFF",
    mint = "#44BB99",
    pear = "#BBCC33",
    olive = "#AAAA00",
    light_yellow = "#EEDD88",
    orange = "#EE8866",
    pink = "#FFAABB",
    pale_grey = "#DDDDDD",
    dark_cyan = "#225555", # added
    dark_blue = "#222255" # added
    )

}

#' States colors 
#'
#' @return a named vector of 8 colors
#' @seealso compute_light_paul_palette()
#' @export
color_states <- function() {

  match_states_colors <- c(
    coexistence = "mint",
    nurse = "pear",
    protegee = "light_cyan",
    desert = "light_yellow",
    protegee_desert = "light_blue",
    nurse_desert = "olive",
    coexistence_desert = "orange",
    protegee_nurse = "pink",
    coexistence_nurse = "dark_cyan",
    coexistence_protegee = "dark_blue",
    unknown = "pale_grey"
    )

  color <- compute_light_paul_palette() %>%
    .[match_states_colors]
  names(color) <- names(match_states_colors)

  return(color)
}

##############
#  Figure 2  #
##############

plot_fig2 <- function(states) {

  #appender <- function(string, suffix = "u = ") { paste0(suffix, string) }
  u_appender <- as_labeller(c("0" = "Without indirect facilitation (0%)", "5" = "High indirect facilitation (40%)", "10" = "Strong indirect facilitation (60%)"))
  stable_states_lab <- c("desert" = "Desert", "protegee_desert" = "Protegee / Desert" , "protegee" = "Protegee", "coexistence" = "Coexistence",
    "coexistence_desert" = "Coexistence / Desert", "nurse_desert" = "Nurse / Desert", "nurse" = "Nurse", "coexistence_nurse" = "Coexistence / Nurse", "coexistence_protegee" = "Coexistence / Protegee")

  g <- plot_diagram(states) +
    xlab(paste("Environmental quality (b)")) +
    ylab(paste("Grazing intensity (g)")) +
    scale_fill_manual(
      labels = as_labeller(stable_states_lab),
      values = color_states(),
      name = "Stable states"
      ) +
    facet_grid(cols = vars(u), labeller = u_appender)# +
    #hrbrthemes::theme_ipsum_rc()

  g
}

plot_bifurcation <- function(scenar, debug_mode = FALSE) {

  #appender <- function(string, suffix = "u = ") { paste0(suffix, string) }
  u_appender <- as_labeller(c("0" = "Without indirect facilitation (0%)", "5" = "High indirect facilitation (40%)", "10" = "Strong indirect facilitation (60%)"))
  stable_states_lab <- c("desert" = "Desert", "protegee_desert" = "Protegee / Desert" , "protegee" = "Protegee", "coexistence" = "Coexistence",
    "coexistence_desert" = "Coexistence / Desert", "nurse_desert" = "Nurse / Desert", "nurse" = "Nurse", "coexistence_nurse" = "Coexistence / Nurse", "coexistence_protegee" = "Coexistence / Protegee")

  # Select variables   
  scenar <- select(scenar, scenario, g, b, u, N, P)
  scenar$run %<>% gather(species, rho, N, P)
  # Make group to split lines
  scenar$run %<>% mutate(group = ifelse(rho > .1, "high", "low")) %>%
    mutate(
      group = as.factor(group),
      species = as.factor(species)
      )

  scenar_high <- filter(scenar, scenario == "together")
  scenar_low <- filter(scenar, scenario == "low_together")

  if (debug_mode) {
   return(scenar_high)
  }

  g <- ggplot2::ggplot(filter(scenar_high$run, group == "high"),
    aes(x = b, y = rho, color = species)) +
    geom_line() +
    ylim(-0.005,1) +
    geom_line(data = filter(scenar_high$run, group == "low"),
      mapping = aes(x = b, y = rho, color = species)) +
    geom_line(data = filter(scenar_low$run, group == "low"),
      mapping = aes(x = b, y = rho - .005, color = species), linetype = "dashed") +
    geom_line(data = filter(scenar_low$run, group == "high"),
      mapping = aes(x = b, y = rho - .005, color = species), linetype = "dashed")
    #xlab(expression(paste("Environmental quality (", b, ")"))) +
    #ylab(expression(paste("Density (", rho, ")"))) #+
    #facet_grid(cols = vars(u), rows = vars(g), labeller = labeller(u = u_appender))# +
    #hrbrthemes::theme_ipsum_rc()

  g
}

plot_fig3 <- function(clustering) {

  g_appender <- function(string, suffix = "g = ") { paste0(suffix, string) }
  clustering$run %<>% gather(var, c, cnn, cnp, cpp, cveg)

  cxx <- as_labeller(c(cnn = "Nurse / Nurse", cnp = "Nurse / Protegee",
      cpp = "Protegee / Protegee", cveg = "Vegetation / Vegetation"))

  g <- plot_diagram(clustering, param = c(x = "del", y = "u"), fill = "c") +
    facet_grid(vars(g), vars(var), labeller = labeller(g = as_labeller(g_appender), var = cxx)) +
    labs(
      x = expression(paste("Proportion of global dispersal (", delta, ")")),
      y = expression(paste("Strength of grazing protection (", u, ")")),
      fill = "Clustering") +
    scale_fill_gradient2(
      #trans = "log10",
      midpoint = 1,
      low = scales::muted("blue"),
      mid = "white",
      high = scales::muted("red"),
      guide = guide_colorbar(title.position = "top")) +
    hrbrthemes::theme_ipsum_rc()

  g
}

plot_fig3bis <- function(clustering, x = "b", y = "g", facet = "u") {

  g_appender <- function(string, suffix = paste(facet, "= ")) { paste0(suffix, string) }
  clustering$run %<>% gather(var, c, cnn, cnp, cpp)#, cveg

  cxx <- as_labeller(c(cnn = "Nurse / Nurse", cnp = "Nurse / Protegee",
      cpp = "Protegee / Protegee", cveg = "Vegetation / Vegetation"))

  g <- plot_diagram(clustering, param = c(x = x, y = y), fill = "c") +
    facet_grid(vars(get(facet)), vars(var), labeller = labeller(g = as_labeller(g_appender), var = cxx)) +
    labs(
      x = paste("Environmental quality  (", x, ")"),
      y = paste("Grazing intensity (", y, ")"),
      fill = "Clustering") +
    scale_fill_gradient2(
      #trans = "log10",
      midpoint = 1,
      low = scales::muted("blue"),
      mid = "white",
      high = scales::muted("red"),
      guide = guide_colorbar(title.position = "top")) +
    hrbrthemes::theme_ipsum_rc()
  g
}

plot_fig3bisbis <- function(clustering, x = "b", y = "g", facet = "u") {

  g_appender <- function(string, suffix = paste(facet, "= ")) { paste0(suffix, string) }
  clustering$run %<>% gather(var, c, cnp)

  cxx <- as_labeller(c(cnn = "Nurse / Nurse", cnp = "Nurse / Protegee",
      cpp = "Protegee / Protegee", cveg = "Vegetation / Vegetation"))

  g <- plot_diagram(clustering, param = c(x = x, y = y), fill = "c") +
    labs(
      x = paste("Fraction of global dispersal (", x, ")"),
      y = paste("Strength of grazing protection (", y, ")"),
      fill = "Clustering") +
    scale_fill_gradient2(
      #trans = "log10",
      midpoint = 1,
      low = scales::muted("blue"),
      mid = "white",
      high = scales::muted("red"),
      guide = guide_colorbar(title.position = "top")) +
    facet_grid(vars(get(facet)), vars(var), labeller = labeller(g = as_labeller(g_appender), var = clustering_labeller())) + theme_alain()
  g
}

theme_alain <- function(){

  hrbrthemes::theme_ipsum_rc() +
  theme(#Whipe
    text = element_text(family = "Helvetica", hjust = .5),
    axis.title = element_text(family = "Helvetica", hjust = .5),
    axis.title.x = NULL,
    axis.title.y = NULL,
    axis.text.x = NULL,
    axis.text.y = NULL,
    strip.text = NULL) +
  theme(#Set up
    axis.title.y = element_text(angle = 90),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 8),
    plot.margin = unit(c(.5, .5, .5, .5), "cm")
    )

}

scale_fill_temperature <- function (mid = 0) {
  scale_fill_gradient2(
      midpoint = mid,
      low = scales::muted("blue"),
      mid = "white",
      high = scales::muted("red"))
}

clustering_labeller <- function () {
  as_labeller(c(cnn = "Nurse / Nurse", cnp = "Nurse / Protegee",
      cpp = "Protegee / Protegee", cveg = "Vegetation / Vegetation"))
}

species_labeller <- function () {
  as_labeller(c(N = "Nurse", P = "Protegee"))
}

facet_labeller <- function (prefix = "b", sep = "= ") {

  prefix <- paste(prefix, sep)
  return(
    as_labeller(function(string, suffix = prefix) { paste0(suffix, string) })   
    
    )
  #g_appender <- function(string, suffix = ) { paste0(suffix, string) }
}

  stable_states_labeller <- function () {
  as_labeller(
    c(
      "desert" = "Desert",
      "protegee_desert" = "Protegee / Desert",
      "protegee" = "Protegee",
      "coexistence" = "Coexistence",
      "coexistence_desert" = "Coexistence / Desert",
      "nurse_desert" = "Nurse / Desert",
      "nurse" = "Nurse",
      "coexistence_nurse" = "Coexistence / Nurse",
      "coexistence_protegee" = "Coexistence / Protegee"
      )
    )
} 

xylabs <- function (...) {
  dots <- pryr::named_dots(...)
  dots <- unlist(dots)

  lab_list <- list(
    del = expression(paste("Fraction of global dispersal (", delta, ")")),
    u = paste("Strength of grazing protection (u)"),
    rho = expression(paste("Species density (", rho, ")")),
    g = paste("Grazing intensity (g)"),
    b = paste("Environmental quality (b)")
    )
  
  lab_used <- lab_list[dots]
  names(lab_used) <- names(dots)

  labs(
    x = lab_used["x"][[1]],
    y = lab_used["y"][[1]]
    )
  #lab_used["x"]
}

