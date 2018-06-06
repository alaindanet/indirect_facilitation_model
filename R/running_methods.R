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
  param = NULL, nb_cores = NULL, solver_type = NULL, inits = NULL) {

  if (is.null(gradientx)) {
    gradientx <- gradienty
  }
  # Prepare the combination of parameters
  gradient <- expand.grid(y = gradienty, x = gradientx) %>%
    tibble::as.tibble(.)
  colnames(gradient) <- c(y, x)
  model <- model_spec
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
	param = simecol::parms(model)[which(!names(simecol::parms(model)) %in% c(x, y))],
	run = output
	),
    class = c("list", "gradient"))
    )
  } else if (identical(run_type, run_bifurc_model, ignore.bytecode = FALSE)){
  return(
    structure(
      list(
	param = simecol::parms(model)[which(!names(simecol::parms(model)) %in% c(x, y))],
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

#' Run simulation over different initial scenarii
#'
#' A wrapper of run_2d_gradient which takes a list of initial values
#'
#' @inheritParams run_2d_gradient
#' @param scenarii a list of vector. A vector contains initial values of the
#' model  
#'
#' @export
run_scenarii_gradient <- function (y = "g", x = "b",
  gradienty = seq(0, .3, by = .1), gradientx = seq(0, 1, by = 0.1),
  model_spec = two_facilitation_model(),
  time_seq = c(from = 0, to = 1000, by = 1),
  param = NULL, nb_cores = NULL, solver_type = NULL,
  scenarii = init_scenarii()){

  # Run the simulations
  output <- tibble::tibble(
    scenario = names(scenarii),
    inits = scenarii) %>%
  dplyr::mutate(gradient = purrr::map(inits, ~ run_2d_gradient(y, x,
	gradienty, gradientx,
	model_spec,
	run_2d_model,
	time_seq,
	param, nb_cores, solver_type, inits = .x))
    )

  # Save model parameters
  model <- model_spec
  basis_param <- simecol::parms(model)[which(!names(simecol::parms(model)) %in%
    c(x, y))]
  if (!is.null(param)) {
    basis_param[names(param)] <- param
  }

  return(
    structure(
      list(
	param = basis_param, run = output),
    class = c("list", "scenarii"))
    )

}
#' TODO: Document function

#' Initialize starting values of state variables
#' 
#' Compute initial values of states variable according to a scenario
#' 
#' @param type character. the scenario
#' @param ini_cover numeric initial cover for high starting cover
#'
#' @details when type = "all", all the scenarii are returned. High cover
#' scenarii set the total cover to ini_cover. It is divided by 2 in together. In
#' low_P or low_N scenarii, the cover of the rare species is of low_cover
#' divided by 2.
#'
#' @return A named list. Each element contains a numeric vector of initial
#' values of the state variables.
#' 
#' @export
init_scenarii <- function (type = "together",
  model = two_facilitation_model(),
  ini_cover = .8, low_cover = .05) {

  stopifnot(type %in% c("nurse", "protegee", "together", "low_N", "low_P", "low_together", "all", "bifurcation"))

  # three or four states model ?
  variables <- names(simecol::init(model))
  if (all(c("N", "P", "D", "ND", "NN", "PP", "NP", "PD", "DD") %in% variables)){
    # four states
    nurse_only <- c(N = ini_cover, P = 0, D = .1,
      NP = 0, PP = 0, NN = ini_cover * ini_cover,
      DD = .1 * .1, PD = 0, ND = ini_cover * .1)
    protegee_only <- c(N = 0, P = ini_cover, D = .1,
      NP = 0, PP = ini_cover * ini_cover, NN = 0,
      DD = .1 * .1, PD = ini_cover * .1, ND = 0)

    mi_cover <- ini_cover / 2
    mi_low_cover <-  low_cover / 2
    high_cover <- ini_cover - mi_low_cover

    low_N <- c(N = mi_low_cover, P = high_cover, D = .1,
      NP = high_cover * mi_low_cover, PP = high_cover * high_cover,
      NN = mi_low_cover * mi_low_cover, DD = .1 * .1,
      PD = high_cover * .1, ND = mi_low_cover * .1)
    low_P <- c(N = high_cover, P = mi_low_cover, D = .1,
      NP = high_cover * mi_low_cover, PP = mi_low_cover * mi_low_cover,
      NN = high_cover * high_cover, DD = .1 * .1,
      PD = mi_low_cover * .1, ND = high_cover * .1)

    low_together <- c(N = mi_low_cover, P = mi_low_cover, D = .1,
      NP = mi_low_cover * mi_low_cover, PP = mi_low_cover * mi_low_cover,
      NN = mi_low_cover * mi_low_cover, DD = .1 * .1,
      PD = mi_low_cover * .1, ND = mi_low_cover * .1)

    together <- c(N = mi_cover, P = mi_cover, D = .1,
      NP = mi_cover * mi_cover, PP = mi_cover * mi_cover,
      NN = mi_cover * mi_cover, DD = .1 * .1,
      PD = mi_cover * .1, ND = mi_cover * .1)


    # TODO: Reframe this: build the list with all init vectors and subset it
    # according to the option provided

    all_inits <- list(
      "protegee_only" = protegee_only,
      "nurse_only" = nurse_only,
      "together" = together,
      "low_N" = low_N,
      "low_P" = low_P,
      "low_together" = low_together
      )
    if(all(type == "all")) {
      return(all_inits)
    } else if (all(type == "bifurcation")) {
      return(all_inits[c("low_together", "together")])
    } else {
      return(all_inits[type])
    }
  } else {
    message("Only the four state model is specified")
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
    dplyr::slice( (n() - cut_row) : n()) %>% # Keep the last simulation
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
avg_runs.gradient <- function(run, cut_row = 10, nb_cores = NULL) {

  param <- run[["param"]]
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
      list(param = param,
      run = run),
    class = c("list", "gradient")
    )
    )
}
avg_runs.scenarii <- function(scenarii, cut_row = 10, nb_cores = NULL) {

  param <- scenarii$param
  run <- scenarii[["run"]] %>%
    dplyr::mutate(
      avg = purrr::map(gradient, avg_runs.gradient, cut_row, nb_cores),
      avg = purrr::map(avg, function(x) x$run)# Keep only the runs
      ) %>%
  unnest(avg)

  return(
    structure(
      list(param = param,
      run = run),
    class = c("list", "scenarii")
    )
    )
}
avg_runs.bifurcation <- function(x, cut_row = 10) {
  output <- avg_runs.gradient(x, cut_row)
  class(output) <- c("bifurcation", "list")
  return(output)
}

#' Convert a gradient object to a scenarii one
#'
#' @param x a gradient object
#' @param scenario a character vector of length one. See @details   
#' @param 
#' @details The scenario should describe one of init density values defined in
#' init_scenarii 
#' @seealso init_scenarii
#' @export
convert2scenarii <- function (x, scenario) {

  basis_param <- x$param

  output <- x$run %>%
    dplyr::mutate(
      scenario = scenario,
      inits = map(scenario, init_scenarii)
      )

  return(
    structure(
      list(
	param = basis_param, run = output),
      class = c("list", "scenarii"))
    )

}

#' Bind scenarii object together 
#' 
#' Combine a list of scenarii in one
#'
#' @param l a list of scenarii object  
#' @param var_name a character vector. Names of the parameters which are different
#' between scenarii
#' @param 
#' @details The scenario should describe one of init density values defined in
#' init_scenarii 
#' @seealso init_scenarii
#' @export
bind_scenar <- function(x, ...) UseMethod("bind_scenar")
bind_scenar.default <- function(x) "Unknown class"
bind_scenar.list <- function (l, var_name = "u") {

  # Parameters
  param <- lapply(l, function(x) x$param)
  common_param <- param[[1]][which(!names(param[[1]]) %in% var_name)]
  varing_param <- sapply(param, function(param)
    param[which(names(param) %in% var_name)])

  if(!all(sapply(varing_param, function(x){ length(x) > 0 })) ){
    warning(paste("There was a problem.", var_name, "was not found to be a varing
	parameter across scenarii."))
  }

  output <- mapply(
  function(param, var_name, scenar) {
    run <- scenar$run
    columns <- lapply(param, function(x) {rep(x, nrow(run))})
    names(columns) <- var_name
    test <- cbind(run, columns)
    return(as.tibble(test))
    },
    varing_param, names(varing_param), l, USE.NAMES = FALSE, SIMPLIFY = FALSE)
  output <- do.call(rbind,
    lapply(output, as.tibble, stringsAsFactors=FALSE)
    )

  return(list(
      param = common_param,
      run = output)
    )

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
  threshold = 10 ^ - 3, multi = FALSE) {

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
def_multi_states <- function(solo, together, sim_status,
  threshold = 10 ^ - 3, multi = FALSE) {
  #TODO: To finish

    if (!sim_status) {
      return("warning")
    } else if (solo == "nurse" && together == "nurse") {
      return("nurse")
    } else if (solo == "nurse" && together == "protégée") {
      return("nurse_protégé")
    } else if (solo == "nurse" && together == "coexistence") {
      return("nurse_coexistence")
    } else if (solo == "nurse" && together == "desert") {
      return("nurse_desert")
    }
}

#' Compute the state result of a range of simulation 
#' 
#' @param data a gradient object
#' @param param
#' @param possibles_states 
#'
#' @export 
compute_states <- function(x, ...) UseMethod("compute_states")
compute_states.gradient <- function (
  data,
  param, 
  possible_states = c("coexistence", "nurse", "protégée", "extinct", "warning")) {

  param <- data$param
  data %<>% .[["run"]]

  var_to_drop <- names(data)[!(names(data) %in% c(param))]

  data %<>%
    dplyr::mutate(
      state = purrr::pmap_chr(
	list(nurse = N, protegee = P, sim_status = status),
	def_state)) %>%
  purrr::modify_at(var_to_drop, ~NULL) %>%
  dplyr::mutate(state = factor(state, levels = possible_states))

return(list(param = param, run = data))

}
compute_states.scenarii <- function (data, param, possible_states = c("coexistence", "nurse", "protégée", "extinct", "warning")) {
  compute_states.gradient(data, param, possible_states)
}
