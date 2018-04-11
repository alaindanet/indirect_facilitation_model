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
  time_seq = c(from = 0, to = 1000, by = 1),
  param = NULL, nb_cores = NULL) {

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
    simecol::parms(model) <- param
  }

  # Run the model
  run_model <- function(x, y) {

    simecol::parms(model)["gamma1"] <- x
    simecol::parms(model)["g"] <- y
    run <- simecol::sim(model)
    output <- simecol::out(run) %>%
      dplyr::select(-time)
    return(output)
  }

  if (is.null(nb_cores)) {
    output <- gradient %>%
      dplyr::mutate(runs = purrr::map2(gamma1, g, run_model))
  } else {
    # In parallel
    cluster <- multidplyr::create_cluster(nb_cores)
    # Export functions of the ODE:
    f_to_load <- c("run_model", "NE_context", "Ncolonize", "PE_context",
      "Pcolonize", "check_nbs", "check_z", "die")
    export_f <- function(x) {
    multidplyr::cluster_assign_value(cluster, x, get(x))
  }
    lapply(f_to_load, export_f)

    multidplyr::set_default_cluster(cluster)
    multidplyr::cluster_library(cluster, c("magrittr"))
    by_var <- gradient %>%
      multidplyr::partition(gamma1) %>%
      dplyr::mutate(runs = purrr::map2(gamma1, g, run_model)) %>%
      dplyr::collect()

    return(by_var)

  }

  return(output)

}

#' Extract the average density of the last timesteps 
#'
#' This function extracts the density of the last timesteps and average it. 
#' 
#' @param run a dataframe. It contains the variable time and the density of the
#' of the state variable.
#' @param n a integer. Define at how many timesteps from the last one should be
#' selected to average density 
#' @return a dataframe. It contains
#' @export
avg_runs <- function(run, cut_row = 100) {

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
  } else if (nurse <= threshold & protegee > threshold) {
    return("protégée")
  } else if (nurse > threshold & protegee <= threshold) {
    return("nurse")
  } else if (nurse <= threshold & protegee <= threshold) {
    return("extinct")
  } else if (nurse > threshold & protegee > threshold) {
    return("coexistence")
  }
}
