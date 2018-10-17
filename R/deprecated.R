################################################################################
#                             Deprecated functions                             #
################################################################################

#####################
#  Running methods  #
#####################

#' Running a model over a two dimension gradient of parameters 
#'
#' This function runs a SimObj model along a two dimensional gradient of
#' parameters   
#' 
#' @param y A string. The first parameter of the 2d gradient
#' @param x A string. The second parameter of the 2d gradient
#' @param gradienty A vector. The vector of parameter values used for the
#' gradient of y. 
#' @param gradientx A vector. Optional. Same as gradienty but for the second
#' parameter. If ommited, the function uses the same gradient than for y.   
#' @param model_spec A character string. This is the name of the model function
#' to use, which calls a SimObj. It can take the values "two_facilitation_model"
#' or "indirect_facilitation_model"
#' @param time_seq A vector. A vector containing the timestep values at which
#' the model will be evaluated.
#' @return A data.frame of tibble data.frame containing the output of the model    
#'
#' @export
run_2d_gradient <- function(y = "g", x = "gamma1",
  gradienty = seq(0, 0.2, length.out = 5), gradientx = NULL,
  model_spec = "indirect_facilitation_model",
  run_type = run_2d_model,
  time_seq = c(from = 0, to = 1000, by = 1),
  param = NULL, nb_cores = NULL, solver_type = NULL, inits = NULL) {

  if (is.null(gradientx)) {
    gradientx <- gradienty
  }
  # Prepare the combination of parameters
  gradient <- expand.grid(y = gradienty, x = gradientx) %>% tibble::as.tibble(.)
  colnames(gradient) <- c(y, x)
  model <- eval(call(model_spec))
  # Define parameters:
  simecol::times(model) <- time_seq
  if (!is.null(param)) {
    simecol::parms(model)[names(param)] <- param
  }
  # Define the solver type
  if (!is.null(solver_type)) {
    simecol::solver(model) <- solver_type
  }
  # Define the initial densities
  if (!is.null(inits)) {
    simecol::init(model) <- inits
  }

  # Run the model
  if (is.null(nb_cores)) {
    output <- gradient %>%
      dplyr::mutate(
	runs = purrr::map2(get(x), get(y), run_type,
	  name_x = x, name_y = y, model = model
	  )
	)
  } else {
    # In parallel
    cluster <- multidplyr::create_cluster(nb_cores)

    # Export functions of the ODE:
    f_to_load <- c("run_type", "model", "NE_context", "Ncolonize", "PE_context",
      "D_context", "Pcolonize", "check_nbs", "check_z", "die",
      "three_states_sys", "four_states_sys", "compute_as", "compute_p_one_z",
      "compute_tau", "degrade", "regen", "facilitate")
    lapply(f_to_load, function(x) {
      multidplyr::cluster_assign_value(cluster, x, get(x))
      })
    # Export arguments of the functions
    multidplyr::cluster_copy(cluster, x)
    multidplyr::cluster_copy(cluster, y)

    multidplyr::set_default_cluster(cluster)
    multidplyr::cluster_library(cluster, c("magrittr", "deSolve"))
    # Run the model
    output <- gradient %>%
      multidplyr::partition() %>%
      dplyr::mutate(
	runs = purrr::map2(
	  get(x),
	  get(y),
	  run_type, name_x = x, name_y = y, model = model)) %>%
      dplyr::collect() %>%
      dplyr::ungroup()
  }
  # Specify output
  if (identical(run_type, run_2d_model, ignore.bytecode = FALSE)) {
  return(
    structure(
      list(
	model = model_spec,
	param = simecol::parms(model)[which(!names(simecol::parms(model)) %in% c(x, y))],
	run = output
	),
    class = c("gradient", "list"))
    )
  } else if (identical(run_type, run_bifurc_model, ignore.bytecode = FALSE)){
  return(
    structure(
      list(
	param = simecol::parms(model)[which(!names(simecol::parms(model)) %in% c(x, y))],
	run = output
	),
    class = c("bifurcation", "list"))
    )
  }
}

#' Run the model by specifying two parameters  
#' 
#' Run the model of a two dimensional gradient.
#' 
#' @param x First variable.
#' @param y Second variable.
#' @param name_x A character vector. 
#' @param name_y A character vector. 
#'
#'
#' @return a data.frame 
#' @export
run_2d_model <- function(x, y, name_x, name_y, model) {

  simecol::parms(model)[name_x] <- x
  simecol::parms(model)[name_y] <- y
  run <- simecol::sim(model)
  output <- simecol::out(run)
  return(output)
}

#' Bifurcation state diagram
#' 
#' This function takes as y the vector of initial species densities and as x the
#' definition of the parameter(s).
#'
#' @details The densities of species pairs are defined as
#' 
#' @inheritParams run_2d_model
#' @seealso init_scenarii
#' @export
run_bifurc_model <- function(x, scenario = "together", name_x, name_y, model) {

  if(any(length(scenario) != 1, !is.character(scenario), scenario == "all")){
    stop("scenario have to be a character vector of length 1")
  }

  simecol::parms(model)[name_x] <- x
  simecol::init(model) <- unlist(
    init_scenarii(type = scenario, model = model)
    )

  run <- simecol::sim(model)
  output <- simecol::out(run) %>%
    dplyr::select(-time)
  return(output)
}

######################
#  Analysis methods  #
######################

avg_runs.gradient <- function(run, cut_row = 10, nb_cores = NULL) {

  param <- run[["param"]]
  model <- run[["model"]]
  run %<>% .[["run"]]

  if (is.null(nb_cores)) {
    run %<>% dplyr::mutate(
      avg = purrr::map(runs, avg_runs, cut_row = cut_row)
      ) %>%
    tidyr::unnest(avg)

  } else {
    cluster <- multidplyr::create_cluster(nb_cores)
    f_to_load <- c("avg_runs.data.frame", "is_run_normal")
    lapply(f_to_load, function(x) {
      multidplyr::cluster_assign_value(cluster, x, get(x))
      })
    # Export arguments of the functions
    multidplyr::cluster_copy(cluster, cut_row)
    multidplyr::cluster_copy(cluster, run)

    multidplyr::set_default_cluster(cluster)
    multidplyr::cluster_library(cluster, c("magrittr", "purrr", "tidyr", "dplyr"))
    # Run the model
    run %<>% multidplyr::partition() %>%
      dplyr::mutate(
	avg = purrr::map(runs, avg_runs.data.frame, cut_row = cut_row),
	avg = purrr::map(avg, clean_run)
	) %>%
    dplyr::collect() %>%
    dplyr::ungroup() %>%
    tidyr::unnest(avg)
  }

  run %<>% clean_run(.)

  return(
    structure(
      list(
	model = model,
	param = param,
	run = run
	),
    class = c("list", "gradient")
    )
    )
}
avg_runs.bifurcation <- function(x, cut_row = 10) {
  output <- avg_runs.gradient(x, cut_row)
  class(output) <- c("bifurcation", "list")
  return(output)
}
compute_occurences.gradient <- function(data, ...) {
  data <- compute_occurences.scenarii(data, ...)

  return(
    structure(
      data,
      class = c("scenarii", "list") #, old_class
      )
    )

}

##################
#  Plot methods  #
##################

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
  ggplot2::geom_raster(interpolate = TRUE) +
  theme_diagram() +
  ggplot2::scale_fill_manual(
    values = col_states,
    limits = possible_states
    )

  return(g)
}
