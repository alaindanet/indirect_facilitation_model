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
      values = color_states(),
      name = "Stable states",
      labels = ggplot2::as_labeller(stable_states_labeller())
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

  if (debug_mode) {
    return(data)
  }
  return(g)
}
plot_diagram.data.frame <- function (
  data,
  param = c(x = "b", y = "g"), debug_mode = FALSE, fill = "N") {

    g <- ggplot2::ggplot(data,
      aes_string(x = param["x"], y = param["y"], fill = fill)) +
    ggplot2::geom_raster() +
    theme_diagram()

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

  stable_states_lab <- c("desert" = "Desert", "protegee_desert" = "Protegee / Desert" , "protegee" = "Protegee", "coexistence" = "Coexistence",
    "coexistence_desert" = "Coexistence / Desert", "nurse_desert" = "Nurse / Desert", "nurse" = "Nurse", "coexistence_nurse" = "Coexistence / Nurse", "coexistence_protegee" = "Coexistence / Protegee")

  g <- plot_diagram(states) +
    xlab(paste("Environmental quality (b)")) +
    ylab(paste("Grazing intensity (g)")) +
    ggplot2::scale_fill_manual(
      labels = ggplot2::as_labeller(stable_states_lab),
      values = color_states(),
      name = "Stable states"
      ) +
    facet_grid(cols = vars(u), labeller = u_labeller)# +
    #hrbrthemes::theme_ipsum_rc()

  g
}

##############
#  FIgure 3  #
##############
identify_discontinuity <- function(x, ...) UseMethod("identify_discontinuity")
identify_discontinuity.default <- function(x, ...) "Unknown class"
identify_discontinuity.numeric <- function(x, threshold = .1) {

#https://stackoverflow.com/a/23863893/5968131
check_diff <- abs(x[-1] - x[-length(x)]) >= threshold
split_at <- function(x, pos) split(x, cumsum(seq_along(x) %in% (pos + 1)))

splitted_list <- split_at(x, which(check_diff))
names(splitted_list) <- 1:length(splitted_list)

l2 <- lapply(seq_along(splitted_list),
  function(y, n, i) {
  as.numeric(rep(n[[i]], length(y[[i]])))
  },
  y = splitted_list, n = names(splitted_list))

as.integer(unlist(l2))
}
identify_discontinuity.data.frame <- function(data, var = "rho", threshold = .1) {

  identify_discontinuity(unlist(data[, var]), threshold = threshold)

}
plot_bifurcation <- function(scenar, ...) {
  sp_var <- rlang::quos(...)

  # Select variables   
  scenar <- dplyr::select(scenar, scenario, g, b, u, !!!sp_var)
  scenar$run %<>%
    dplyr::mutate(P = P - 0.005) %>%
    tidyr::gather(species, rho, !!!sp_var) %>%
    dplyr::mutate(species = as.factor(species))
  # Make group to split lines
  scenar$run %<>%
    dplyr::group_by(scenario, g, u, species) %>%
    dplyr::group_nest() %>%
    dplyr::mutate(
      group = purrr::map(data, identify_discontinuity, var = "rho", threshold = .05)
    ) %>%
    tidyr::unnest()

  # Have to do this in two steps due to a ggplot2 limitation: cannot specify
  # group and linetype in the same time:
  #https://stackoverflow.com/a/27011361/5968131
  ## Create the two:
  scenario_list <- lapply(names(scenar$inits), function(x) {
    filter(scenar$run,
      scenario == x)
      })
  g <- ggplot2::ggplot(scenario_list[[2]],
    aes(x = b, y = rho, color = species,
      group = interaction(group, species))) +
    geom_line() +
    ylim(-0.005, 1)

  g + geom_line(mapping = aes(y = rho, group = interaction(group, species)),#y = rho - 0.005
    data = scenario_list[[1]], linetype = "solid")
}

theme_alain <- function(){

  hrbrthemes::theme_ipsum_rc() +
  theme(#Whipe
    text = element_text(family = "Helvetica", hjust = .5),
    axis.title = element_text(family = "Helvetica", hjust = .5, face = "bold", size = 10),
    plot.title = NULL,
    axis.title.x = NULL,
    axis.title.y = NULL,
    axis.text.x = NULL,
    axis.text.y = NULL,
    strip.text = NULL) +
  theme(#Set up
    plot.title = element_text(family = "Helvetica", hjust = .5, face = "bold", size = 10),
    axis.title.y = element_text(angle = 90, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 8, margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")),
    plot.margin = unit(c(.5, .5, .5, .5), "pt"),
    panel.spacing = unit(5, "points")
    )

}

scale_fill_temperature <- function (colors = c("white", "yellow", "orange", "red"), name = NULL, limits = c(NA, NA)) {
  scale_fill_gradientn(
    colors = colors,
    name = name,
    limits = limits
    )
}

scale_colour_species <- function(){
  ggplot2::scale_colour_manual(
    name = "Species",
    values = c(
      N = "#BBCC33",
      P = "#99DDFF",
      total = "black"),
    labels = species_labeller()
    )
}

scale_linetype_model <- function(){
  ggplot2::scale_linetype_manual(
    name = "Model",
    values = c(ca = "solid", pa = "dashed"),
    labels = model_labeller()
    )
}

#######################################################################################################
#  Create own palette:
#  https://drsimonj.svbtle.com/creating-corporate-colour-palettes-for-ggplot2  #
#######################################################################################################

temp_colors <- c(
  `red`        = "#d11141",
  `green`      = "#00b159",
  `blue`       = "#00aedb",
  `orange`     = "#f37735",
  `yellow`     = "#ffc425",
  `light grey` = "#cccccc",
  `dark grey`  = "#8c8c8c")

temp_cols <- function(...) {
  cols <- c(...)

  if (is.null(cols))
    return (temp_colors)

  temp_colors[cols]
}

temp_palettes <- list(
  `main`  = temp_cols("blue", "green", "yellow"),

  `cool`  = temp_cols("blue", "green"),

  `hot`   = temp_cols("yellow", "orange", "red"),

  `mixed` = temp_cols("blue", "green", "yellow", "orange", "red"),

  `grey`  = temp_cols("light grey", "dark grey")
)

temp_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- temp_palettes[[palette]]

  if (reverse) pal <- rev(pal)

  colorRampPalette(pal, ...)
}

## Auto lab:
clustering_labeller <- function () {
  ggplot2::as_labeller(c(cnn = "Nurse / Nurse", cnp = "Nurse / Protegee",
      cpp = "Protegee / Protegee", cveg = "Vegetation / Vegetation"))
}

species_labeller <- function () {
  ggplot2::as_labeller(
    c(N = "Nurse",
      P = "Protegee",
      total = "Total")
    )
}

model_labeller <- function () {
  ggplot2::as_labeller(
    c(ca = "Cellular automata", pa = "Pair approximation")
    )
}

facet_labeller <- function (prefix = "b", sep = "= ") {

  prefix <- paste(prefix, sep)
  return(
    ggplot2::as_labeller(function(string, suffix = prefix) { paste0(suffix, string) })   
    
    )
  #g_appender <- function(string, suffix = ) { paste0(suffix, string) }
}



u_labeller <- ggplot2::as_labeller(c(
    "0" = "Without indirect facilitation (u = 0)",
    "5" = "High indirect facilitation (u = 5)",
    "10" = "Strong indirect facilitation (u = 10)"))

stable_states_labeller <- function () {
  ggplot2::as_labeller(
    c(
      "desert" = "Desert",
      "protegee_desert" = "Protegee / Desert",
      "protegee" = "Protegee",
      "coexistence" = "Coexistence",
      "coexistence_desert" = "Coexistence / Desert",
      "nurse_desert" = "Nurse / Desert",
      "nurse" = "Nurse",
      "coexistence_nurse" = "Coexistence / Nurse",
      "coexistence_protegee" = "Coexistence / Protegee",
      "unknown" = "Inconsistent"
      )
    )
} 

xylabs <- function (...) {
  dots <- pryr::named_dots(...)
  dots <- map(dots, eval) %>% unlist

  lab_list <- list(
    del = expression(bold(paste("Fraction of global dispersal (", delta, ")"))),
    u = paste("Strength of grazing protection (u)"),
    rho = expression(bold(paste("Species density (", rho, ")"))),
    g = paste("Grazing intensity (g)"),
    b = paste("Environmental quality (b)"),
    f = paste("Strength of direct facilitation (f)"),
    time = paste("Time step")
    )
  
  lab_used <- lab_list[dots]
  names(lab_used) <- names(dots)

  labs(
    x = lab_used["x"][[1]],
    y = lab_used["y"][[1]]
    )
}

##############
#  Figure 3  #
##############

plot_fig3 <- function(clustering) {

  g_appender <- function(string, suffix = "g = ") { paste0(suffix, string) }
  clustering$run %<>% gather(var, c, cnn, cnp, cpp, cveg)

  cxx <- ggplot2::as_labeller(c(cnn = "Nurse / Nurse", cnp = "Nurse / Protegee",
      cpp = "Protegee / Protegee", cveg = "Vegetation / Vegetation"))

  g <- plot_diagram(clustering, param = c(x = "del", y = "u"), fill = "c") +
    facet_grid(vars(g), vars(var), labeller = labeller(g = ggplot2::as_labeller(g_appender), var = cxx)) +
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
  clustering$run %<>% tidyr::gather(var, c, cnn, cnp, cpp)#, cveg

  cxx <- ggplot2::as_labeller(c(cnn = "Nurse / Nurse", cnp = "Nurse / Protegee",
      cpp = "Protegee / Protegee", cveg = "Vegetation / Vegetation"))

  g <- plot_diagram(clustering, param = c(x = x, y = y), fill = "c") +
    facet_grid(vars(get(facet)), vars(var), labeller = labeller(g = ggplot2::as_labeller(g_appender), var = cxx)) +
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

  cxx <- ggplot2::as_labeller(c(cnn = "Nurse / Nurse", cnp = "Nurse / Protegee",
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
    facet_grid(vars(get(facet)), vars(var), labeller = labeller(g = ggplot2::as_labeller(g_appender), var = clustering_labeller())) + theme_alain()
  g
}

##########################
#  Clustering + contour  #
##########################
plot_clustering_contour <- function (clust, fac = u) {
  fac <- rlang::enquo(fac)

  # Plot by clust_type
  clust %<>%
    dplyr::group_by(clust_type) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      plot_diagram = purrr::map(data,
	plot_clustering, yvar = rlang::quo_name(fac))
    )
      brks_clustering <- c(3, 3.5, 4, 4.1, 4.2, 4.5)
      clust %<>%
	dplyr::mutate(interp = purrr::map2(data, plot_diagram, function(data, graph, brks){
	    #Interp
	    interp <- interp_contour(data, x = del, y = !!fac, z = clustering,
	      nb_pts = 30, duplic = clust_type)
	    # Draw contours
	    contours <- graph +
	      draw_contour(interp, x = x, y = y, z = clustering,
		colour = "black", breaks = brks_clustering)
	    direct.label(contours, list("bottom.pieces", colour = "black"))

}, brks = brks_clustering))
      clust
}

#' Plot clustering 
#'
#' @param clust_data clustering dataset
#' @param yvar character
#'
#'
plot_clustering <- function(clust_data, yvar) {
  ggplot2::ggplot(clust_data,
    aes_string(
      x = "del", y = yvar, fill = "clustering")) +
  ggplot2::geom_raster() +
  theme_alain() +
  scale_fill_temperature(name = "Clustering", limits = c(2, 4.5)) 
}


#' interp contours
#'
#'
interp_contour <- function(clust, ...) UseMethod("interp_contour")
interp_contour.default <- function(clust, ...) "Unknown class"
interp_contour.data.frame <- function (clust, x = NULL, y = NULL, z = NULL, nb_pts = 100, duplic = NULL, ...) {
  var_fill <- rlang::enquo(z)
  y <- rlang::enquo(y)
  x <- rlang::enquo(x)
  duplic <- rlang::enquo(duplic)
  #plot_arg <- rlang::quos(...)

  #Interp
  clust_interp <- na.omit(clust)
  # Check if duplicated row
  if (dplyr::distinct(dplyr::select(clust_interp, !!x, !!y)) %>% nrow /
    nrow(clust_interp) < 1) {

    stopifnot(!rlang::quo_is_null(duplic))
    type <- unique(clust_interp[, rlang::quo_name(duplic)]) %>% unlist 

    for (i in 1:length(type)) {
      temp <- dplyr::filter(clust_interp, !!duplic == type[i])

      interp_temp <- akima::interp(
	x = temp[, rlang::quo_name(x)] %>% unlist,
	y = temp[, rlang::quo_name(y)] %>% unlist,
	z = temp[, rlang::quo_name(var_fill)] %>% unlist,
	xo = seq(
	  min(temp[, rlang::quo_name(x)] %>% unlist),
	  max(temp[, rlang::quo_name(x)] %>% unlist),
	  length.out = nb_pts),
	yo = seq(
	  min(temp[, rlang::quo_name(y)] %>% unlist),
	  max(temp[, rlang::quo_name(y)] %>% unlist),
	  length.out = nb_pts),
	#jitter.random = TRUE,
	linear = TRUE
      )

      interp_temp_df <- akima::interp2xyz(interp_temp, data.frame = TRUE) %>%
	as.tibble %>%
	rename(!!var_fill := z) %>%
	mutate(!!duplic := type[i])
      if (i == 1) {
	interp_df <- interp_temp_df
      } else {
	interp_df <- rbind(interp_df, interp_temp_df)
      }
    }
  } else {
    interp_akima <- akima::interp(
      x = clust_interp[, rlang::quo_name(x)] %>% unlist,
      y = clust_interp[, rlang::quo_name(y)] %>% unlist,
      z = clust_interp[, rlang::quo_name(var_fill)] %>% unlist,
      xo = seq(
	min(clust_interp[, rlang::quo_name(x)] %>% unlist),
	max(clust_interp[, rlang::quo_name(x)] %>% unlist),
	length.out = nb_pts),
      yo = seq(
	min(clust_interp[, rlang::quo_name(y)] %>% unlist),
	max(clust_interp[, rlang::quo_name(y)] %>% unlist),
	length.out = nb_pts),
      #jitter.random = TRUE,
      linear = TRUE
    )
    interp_df <- akima::interp2xyz(interp_akima, data.frame = TRUE) %>%
      as.tibble %>%
      rename(!!var_fill := z)
  }
  interp_df
  # Contour
}

interp_contour.avg_scenarii <- function (clust, x = NULL, y = NULL, z = NULL, nb_pts = 100, duplic = NULL, ...) {
  interp_contour.data.frame(clust$run, x, y, z, nb_pts, duplic, ...)
}

#' draw contours
#'
#'
draw_contour <- function (data, x = NULL, y = NULL, z = NULL, ...) {
  var_fill <- rlang::enquo(z)
  y <- rlang::enquo(y)
  x <- rlang::enquo(x)

  ggplot2::stat_contour(data = data,
    aes(x = !!x, y = !!y, z = !!var_fill, colour = ..level..), ...)
  #Strange error with directlabels: have to specify colour = ..level..
  #https://rdrr.io/rforge/directlabels/src/etc/contour.R
}
#' Remove x axis
#'
#' @param x a ggplot object
rm_x_axis <- function (x) {
  x + theme(axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())
}
rm_y_axis <- function (x) {
  x + theme(axis.title.y = element_blank())
}
