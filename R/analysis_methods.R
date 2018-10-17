
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

  if (nrow(run) < cut_row) {
    warnings("cut_row is ", cut_row, "  and the number of row of the data.frame
      is ", nrow(run), ". cut_row has been set to ", nrow(run), ".")
    cut_row <- nrow(run)
  }

  out <- run %>%
    dplyr::slice( (n() - cut_row) : n()) %>% # Keep the last simulation
    tidyr::gather(species, rho, -time) %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(
      rho = mean(rho, na.rm = TRUE),
      time = last(time)
      ) %>%
    dplyr::ungroup() %>%
    tidyr::spread(species, rho)

  # check negative or NaN
  out$status <- is_run_normal(run)

  return(out)
}
avg_runs.scenarii <- function(scenarii, cut_row = 10) {

  run <- scenarii[["run"]] %>%
    dplyr::mutate(
      avg = parallel::mclapply(run, avg_runs, cut_row = cut_row)
      ) %>%
  tidyr::unnest(avg) %>%
  dplyr::select(-run)

  return(
    structure(
      list(
	model = scenarii$model,
	inits = scenarii$inits,
	param = scenarii$param,
	gradient = scenarii$gradient,
	run = run
      ),
    class = c("avg_scenarii","scenarii", "list")
    )
    )
}

#' Compute the co-occurences between nurse and protegee 
#'
#' This function computes the density of cooccuring species.
#' It should be after avg_runs
#' 
#' @param data a dataframe. 
#' @return a dataframe.
#' @export
compute_occurences <- function(x, ...) UseMethod("compute_occurences")
compute_occurences.default <- function(x) "Unknown class"
compute_occurences.avg_scenarii <- function(data, ...) {

  if (data$model != "ca_two_facilitation_model") {
    run <- data[["run"]] %>%
      dplyr::mutate(#qj|i = pij / pi; cij = qj|i / pj 
      cnp = NP / (N * P),
      cnn = NN / (N * N),
      cpp = PP / (P * P),
      cveg = (NN + 2 * NP + PP) / ( (N + P)^2 )

      )
  message("c_veg has a bad formula")
  } else {
    run <- data[["run"]] %>%
      dplyr::mutate(
	cnp = qnp / N,
	cpn = qpn / P,
	cnn = qnn / N,
	cpp = qpp / P,
	cveg = qveg / (N + P)
	)
  }

  return(
    structure(
      list(
	model = data$model,
	inits = data$inits,
	param = data$param,
	gradient = data$gradient,
	run = run
	),
      class = c("avg_scenarii", "scenarii", "list")
      )
    )
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
      inits = purrr::map(scenario, init_scenarii)
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
#' @param ... scenarii objects  
#' @details The scenario should describe one of init density values defined in
#' init_scenarii 
#' @seealso init_scenarii
#' @export
bind_scenar <- function(x, ...) UseMethod("bind_scenar")
bind_scenar.default <- function(x) "Unknown class"
bind_scenar.scenarii <- function (...) {

  scenar <- list(...)

  # Parameters

  param <- lapply(scenar, function(x) x$param)
  test_param <- sapply(param, FUN = identical, param[[1]])
  if (!all(test_param)) {
    warning(
      paste("There was a problem. \n
	Parameters were not found to be the same across scenarii. \n
	Only gradient is expected to vary \n")
	)
  }

  inits <- lapply(scenar, function(x) x$inits)
  if (!all(sapply(inits, FUN = identical, inits[[1]]))) {
    warning(
      paste("There was a problem. inits were not found to be the same across scenarii. \n
	Only gradient is expected to vary")
	)
  }

  binded_run <- dplyr::bind_rows(lapply(scenar, function(x) x$run))

  binded_scenar <- scenar[[1]]
  rm(scenar)
  binded_scenar$run <- binded_run
  rm(binded_run)
  #We need to update the gradient slot
  binded_scenar[["gradient"]] <- binded_scenar[["run"]] %>%
    .[, names(.) %in% names(binded_scenar[["gradient"]])] %>%
    as.list(.) %>%
    lapply(., function(x) unique(x))

  if (any(class(binded_scenar) %in% "avg_scenarii")) {
    class_returned <- c("avg_scenarii","scenarii", "list")
  } else {
    class_returned <- c("scenarii", "list")
  }

  return(
    structure(
      binded_scenar,
      class = class_returned
      )
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
  threshold = 10 ^ - 3) {

  if (!sim_status) {
    return("warning")
  } else if (nurse <= threshold && protegee > threshold) {
    return("protegee")
  } else if (nurse > threshold && protegee <= threshold) {
    return("nurse")
  } else if (nurse <= threshold && protegee <= threshold) {
    return("desert")
  } else if (nurse > threshold && protegee > threshold) {
    return("coexistence")
  }

}
def_multi_state <- function(scenar_one, scenar_two, sim_status,
  possible_states = c("coexistence", "nurse", "protégée", "extinct", "warning"),
  threshold = 10 ^ - 3) {

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
#' @param param parameters of the gradient
#' @param var_to_keep variable name that you want to keep, e.g. scenario for
#' scenarii object (set by default) 
#' @param type character string. The type of state to be computed, either
#' "single" or "double" for see bistable areas 
#' @param warning_as_desert TRUE by default. Convert anormal simulation to
#' desert (convenient for plotting)
#' @seealso def_state
#'
#' @export 
compute_states <- function(x, ...) UseMethod("compute_states")
compute_states.avg_scenarii <- function (
  data,
  var_to_keep = "scenario",
  type = "single", warning_as_desert = TRUE) {

  stopifnot(type %in% c("single", "double"))

  common_param <- data$param
  gradient_names <- names(data$gradient)
  model <- data$model
  run <- data[["run"]]

  var_to_drop <- names(run)[!(names(run) %in% c(gradient_names, var_to_keep))]

  run %<>%
    dplyr::mutate(
      state = purrr::pmap_chr(
	list(nurse = N, protegee = P, sim_status = status),
	def_state)) %>%
  purrr::modify_at(var_to_drop, ~NULL) %>%
    dplyr::mutate(state = as.factor(state))

  if (warning_as_desert) {
    run %<>%
      dplyr::mutate(
	state = stringr::str_replace(state, "warning", "desert"),
	state = as.factor(state)
	)
    message("warning state has been silently replaced by desert")

  }

  if (type == "single") {

  } else if (type == "double") {

    #TODO: to generalize to any scenario
    run %<>%
      dplyr::mutate(state = as.character(state)) %>%
      tidyr::spread(scenario, state) %>%
      dplyr::mutate(
	state = purrr::map2_chr(low_together, together, define_double_state)
	) %>%
      dplyr::select(-low_together, -together)
      #tidyr::unite(state, -param)

  }

  return(
    structure(
      list(
	model = data$model,
	inits = data$inits,
	param = data$param,
	gradient = data$gradient,
	run = run
      ),
    class = c("states_scenarii","scenarii", "list")
    )
    )
}
compute_states.scenarii <- function (data, param, var_to_keep = "scenario",
  type = "single") {

  compute_states.gradient(
    data,
    param = c(param),
    var_to_keep = var_to_keep,
    type = type)

}

#' Define bistable states 
#'
#' @param scenar1 a character string. It is the name the first scenario   
#' @param scenar2 a character string. It is the name the second scenario   
#'
#' @seealso
#' @export
define_double_state <- function (scenar1, scenar2) {

  if (scenar1 == scenar2){
    return(scenar1)
  } else if (all(scenar1 %in% c("protegee", "desert"), scenar2 %in% c("protegee", "desert"))) {
    return("protegee_desert")
  } else if (all(scenar1 %in% c("nurse", "desert"), scenar2 %in% c("nurse", "desert"))) {
    return("nurse_desert")
  } else if (all(scenar1 %in% c("coexistence", "desert"), scenar2 %in% c("coexistence", "desert"))) {
    return("coexistence_desert")
  } else if (all(scenar1 %in% c("protegee", "nurse"), scenar2 %in% c("protegee", "nurse"))) {
    return("protegee_nurse")
  } else if (all(scenar1 %in% c("coexistence", "nurse"), scenar2 %in% c("coexistence", "nurse"))) {
    return("coexistence_nurse")
  } else if (all(scenar1 %in% c("coexistence", "protegee"), scenar2 %in% c("coexistence", "protegee"))) {
    return("coexistence_protegee")
  } else {
    return("unknown")
  }

}

#' Subset runs
#'
#' @param data a scenarii object
#' @param ... filtering rules
#' @details https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html
#' @return a scenarii object
#' @export
filter.scenarii <- function (data, ...) {
  run <- data[["run"]]
  filter_set <- rlang::quos(...)

   run %<>% dplyr::filter(., !!! filter_set)

  data[["run"]] <- run
  # We need to update the gradient slot
  data[["gradient"]] <- data[["run"]] %>%
    .[, names(.) %in% names(data[["gradient"]])] %>%
    as.list(.) %>%
    lapply(., function(x) unique(x))


  if (any(class(data) %in% "avg_scenarii")) {
    class_returned <- c("avg_scenarii","scenarii", "list")
  } else {
    class_returned <- c("scenarii", "list")
  }

  return(
    structure(
      data,
      class = class_returned
      )
    )
}

#' Select variables of runs
#'
#' @param data a scenarii object
#' @param ... variables selected
#' @details https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html
#' @return a scenarii object
#' @seealso dplyr::select 
#' @export
select.scenarii <- function (data, ...) {
  run <- data[["run"]]
  set <- rlang::quos(...)

  param_gradient <- names(data[["gradient"]])
  run %<>% dplyr::select(., scenario, param_gradient,!!! set)

  data[["run"]] <- run
  data[["gradient"]] <- data[["run"]] %>%
    .[, names(.) %in% names(data[["gradient"]])] %>%
    as.list(.)

  if (any(class(data) %in% "avg_scenarii")) {
    class_returned <- c("avg_scenarii","scenarii", "list")
  } else {
    class_returned <- c("scenarii", "list")
  }

  return(
    structure(
      data,
      class = class_returned
      )
    )
}

#' Check consistency of cellular automata runs 
#'
#' @param data a data.frame 
#' @details It checks that more than 80% of the simulations give the same
#' behaviour (e.g. same species win)
#' @return logical 
#' @seealso
#' @export
check_consistency <- function(data, threshold = .001) {
  # possible outputs
  outcome <- list(
  N_win = length(which(data$N > threshold & data$P < threshold)),
  P_win = length(which(data$P > threshold & data$N < threshold)),
  N_P_win = length(which(data$P > threshold & data$N > threshold)),
  no_win = length(which(data$P < threshold & data$N < threshold))
  )
  # nb replicates:
  nb_rep <- nrow(data)

  test <- sapply(outcome, function(x) if(x/nb_rep >= .8){TRUE}else{FALSE})
  # TODO: return true or false for each case, compute the specific clustering
  #only if the corresponding variables is true

  if(any(test)) {
    return(TRUE)
  } else {
    return(FALSE)
  }

}

#' Identify abrupt transitions
#'
#' @param v a vector 
#' @param threshold numeric 
#' @details It identifies changes in variable values and assign series
#' accordingly to groups. 
#' @return factor 
#' @seealso plot_bifurcation
#' @export
identify_transition <- function(data, threshold = .001) {

  # Which are different ?
  ind <- abs(
  as.numeric(v[-1]) - as.numeric(v[-length(v)])
  ) >= threshold

  # Assign to group
  splitAt <- function(x, pos) split(x, cumsum(seq_along(x) %in% (pos+1)))
  l1 <- splitAt(as.numeric(df_high$run$N), which(ind))
  names(l1) <- 1:length(l1)
  l2 <- lapply(seq_along(l1), function(y, n, i) {
	       as.numeric(rep(n[[i]], length(y[[i]])))
                               }, 
	       y = l1,
	       n = names(l1)
	       )
  unlist(l2)
  #TODO: Test the function

}
