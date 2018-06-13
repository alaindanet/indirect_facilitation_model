
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
    return("protégée")
  } else if (nurse > threshold && protegee <= threshold) {
    return("nurse")
  } else if (nurse <= threshold && protegee <= threshold) {
    return("extinct")
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
#' @param param
#' @param type   
#' @param possibles_states 
#'
#' @export 
compute_states <- function(x, ...) UseMethod("compute_states")
compute_states.gradient <- function (
  data,
  param,
  var_to_keep = NULL,
  possible_states = c("coexistence", "nurse", "protégée", "extinct", "warning"),
  type = "single") {

  common_param <- data$param
  data %<>% .[["run"]]

  var_to_drop <- names(data)[!(names(data) %in% c(param, var_to_keep))]


    data %<>%
      dplyr::mutate(
	state = purrr::pmap_chr(
	  list(nurse = N, protegee = P, sim_status = status),
	  def_state)) %>%
    purrr::modify_at(var_to_drop, ~NULL) %>%
    dplyr::mutate(state = factor(state, levels = possible_states))

  if (type == "single") {


  } else if (type == "double") {

    data %<>%
      tidyr::spread(scenario, state) %>%
      tidyr::unite(states, -param)

  }

return(
  list(param = common_param, run = data)
  )
}

compute_states.scenarii <- function (data, param, var_to_keep = "scenario",
  possible_states = c("coexistence", "nurse", "protégée", "extinct", "warning"),
  type = "single") {

  compute_states.gradient(
    data,
    param = c(param),
    var_to_keep = var_to_keep,
    possible_states,
    type = type)

}
