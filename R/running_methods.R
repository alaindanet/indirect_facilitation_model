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
  time_seq = c(from = 0, to = 100, by = .1)) {

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

  # Run the model
  run_model <- function(x, y) {

    simecol::parms(model)["gamma1"] <- x
    simecol::parms(model)["g"] <- y
    run <- simecol::sim(model)
    output <- simecol::out(run) %>%
      dplyr::select(-time)
    return(output)

  }
  # Select the last values 

  # TODO: generalize for any argument
  output <- gradient %>%
    dplyr::mutate(runs = purrr::map2(gamma1, g, run_model))
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
  # Keep the last 100
  out  <- run %>%
    dplyr::slice(n() - cut_row: n()) %>%
    tidyr::gather(species, rho) %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(rho = mean(rho)) %>%
    dplyr::ungroup() %>%
    tidyr::spread(species, rho)

  return(out)
}
