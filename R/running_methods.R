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
#' @param model_spec A SimObj. A model of the object class defined in the
#' simecol package.
#' @param time_seq A vector. A vector containing the timestep values at which
#' the model will be evaluated.
#' @return A data.frame of tibble data.frame containing the output of the mod#' el    
#'
#' @export
run_2d_gradient <- function(y = "g", x = "gamma1",
  gradienty = seq(0, 0.2, length.out = 5), gradientx = NULL,
  model_spec = indirect_facilitation_model(),
  run_type = run_2d_model,
  time_seq = c(from = 0, to = 1000, by = 1),
  param = NULL, nb_cores = NULL, solver_type = NULL) {

  if (is.null(gradientx)) {
    gradientx <- gradienty
  }
  # Prepare the combination of parameters
  gradient <- expand.grid(y = gradienty, x = gradientx) %>%
    tibble::as.tibble(.)
  colnames(gradient) <- c(y, x)
  # Prepare the model:
  model <- model_spec
  simecol::times(model) <- time_seq
  if (!is.null(param)) {
    simecol::parms(model)[names(param)] <- param
  }
  # Define the solver type
  if (!is.null(solver_type)) {
    simecol::solver(model) <- solver_type
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
	param = parms(model)[which(!names(parms(model)) %in% c(x, y))],
	run = output
	),
    class = c("list", "gradient"))
    )
  } else if (identical(run_type, run_bifurc_model, ignore.bytecode = FALSE)){
  return(
    structure(
      list(
	param = parms(model)[which(!names(parms(model)) %in% c(x, y))],
	run = output
	),
    class = c("list", "bifurcation"))
    )
  }
}

#' Run bifurcation over a gradient
#'
#' A wrapper of run_2d_gradient 
#'
#' @inheritParams run_2d_gradient
#'
#' @export
run_bifurcation <- function(y = "init", x = "b", gradienty = c(.05, .4),
  gradientx = seq(0, 1, by = 0.1), model_spec = indirect_facilitation_model(),
  time_seq = c(from = 0, to = 1000, by = 1), param = NULL, nb_cores = NULL,
  solver_type = NULL) {

  run_2d_gradient(y, x,
  gradienty, gradientx,
  model_spec,
  run_bifurc_model,
  time_seq,
  param, nb_cores, solver_type)

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
  output <- simecol::out(run) %>%
    dplyr::select(-time)
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
#' @export
run_bifurc_model <- function(x, y, name_x, name_y, model) {

  simecol::parms(model)[name_x] <- x
  simecol::init(model)[c("N", "P")] <- y

  pair_names <- ! names(simecol::init(model)) %in% c("N", "P")

  simecol::init(model)[pair_names] <- y * .25
  run <- simecol::sim(model)
  output <- simecol::out(run) %>%
    dplyr::select(-time)
  return(output)
}

#' Extract the average density of the last timesteps 
#'
#' This function extracts the density of the last timesteps and average it. 
#' 
#' @param run a dataframe. It contains the variable time and the density of the
#' of the state variable.
#' @param cut_row a integer. Define at how many timesteps from the last one should be
#' selected to average density 
#' @return a dataframe.
#' @export
avg_runs <- function(x, ...) UseMethod("avg_runs")
avg_runs.default <- function(x) "Unknown class"
avg_runs.data.frame <- function(run, cut_row = 10) {

  if (nrow(run) <= cut_row) {
    # Create a NA data.frame of length one
    out <- lapply(vector("list", ncol(run)), function(x) return(NA))
    names(out) <- names(run)
    out %<>% as.data.frame(.)

  } else {

  out <- run %>%
    dplyr::slice( (n() - cut_row) : n()) %>% # Keep the last 100
    tidyr::gather(species, rho) %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(rho = mean(rho, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    tidyr::spread(species, rho)
  }
  # check negative or NaN
  out$status <- is_run_normal(run)

  return(out)
}
avg_runs.gradient <- function(run, cut_row = 100) {

  param <- run[["param"]]
  run %<>% .[["run"]]
  run %<>% dplyr::mutate(
    avg = purrr::map(runs, avg_runs, cut_row = cut_row)) %>%
    tidyr::unnest(avg)

  return(
    structure(
      list(param = param,
      run = run),
    class = c("list", "gradient")
    )
    )
}
avg_runs.bifurcation <- function(x, cut_row = 10) {
  output <- avg_runs.gradient(x, cut_row)
  class(output) <- c("bifurcation", "list")
  return(output)
}

#' Define the state result of a simulation 
#'
#' This function gives the result of a simulation, i.e., "warning", "nurse",
#' "protegee", "coexistence"
#' 
#' @param nurse a dbl of length 1. The density of the nurse at the steady
#' state. 
#' @param protegee a dbl of length 1. The density of the protegee at the steady
#' state. 
#' @return a chr vector of length 1
#' @export
def_state <- function (nurse, protegee, sim_status,
  threshold = 10 ^ - 3) {

  if (!sim_status) {
    return("warning")
  } else if (nurse <= threshold && protegee > threshold) {
    return("protégée")
  } else if (nurse > threshold && protegee <= threshold) {
    return("nurse")
  } else if (nurse <= threshold && protegee <= threshold) {
    return("extinct")
  } else if (nurse > threshold && protegee > threshold) {
    return("coexistence")
  }
}
