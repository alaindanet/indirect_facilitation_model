
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
compute_occurences.scenarii <- function(data, ...) {

  common_param <- data$param
  model <- data$model
  old_class <- class(data)
  data %<>% .[["run"]] %>%
    dplyr::mutate(cnp = NP / N)

  return(
    structure(
    list( model = model,
      param = common_param,
      run = data
      ),
    class = c("scenarii", "list") #, old_class
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

  return(
    structure(
      list(
	param = common_param,
	run = output),
    class = c("list", "scenarii")
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
compute_states.gradient <- function (
  data,
  param,
  var_to_keep = NULL,
  type = "single", warning_as_desert = TRUE) {

  common_param <- data$param
  model <- data$model
  data %<>% .[["run"]]

  var_to_drop <- names(data)[!(names(data) %in% c(param, var_to_keep))]

  data %<>%
    dplyr::mutate(
      state = purrr::pmap_chr(
	list(nurse = N, protegee = P, sim_status = status),
	def_state)) %>%
  purrr::modify_at(var_to_drop, ~NULL) %>%
    dplyr::mutate(state = as.factor(state))

  if (warning_as_desert) {
    data %<>%
      dplyr::mutate(
	state = stringr::str_replace(state, "warning", "desert"),
	state = as.factor(state)
	)
    message("warning state has been silently replaced by desert")


  }

  if (type == "single") {

  } else if (type == "double") {

    #TODO: to generalize to any scenario
    data %<>%
      dplyr::mutate( state = as.character(state)) %>%
      tidyr::spread(scenario, state) %>%
      dplyr::mutate(
	state = purrr::map2_chr(low_together, together, define_double_state)
	) %>%
      dplyr::select(-low_together, -together)
      #tidyr::unite(state, -param)

  }

return(
  list(
    model = model,
    param = common_param,
    run = data
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
#' @seealso a scenarii object
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

