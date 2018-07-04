
#' Run simulation over different initial scenarii
#'
#' A wrapper of run_2d_gradient which takes a list of initial values
#'
#' @inheritParams run_2d_gradient
#' @param scenarii a list of vector. A vector contains initial values of the
#' model  
#'
#' @export
run_scenarii_gradient <- function (
  gradient = list(g = seq(0, .3, length.out = 1),
    b = seq(0, 1, length.out = 1)),
  model_spec = "two_facilitation_model",
  time_seq = c(from = 0, to = 1000, by = 1),
  param = NULL, nb_cores = NULL, solver_type = NULL,
  scenarii = NULL) {

  if (is.null(gradient)) {
    stop("Please provide a gradient")
  }

  model <- eval(call(model_spec))
  # Define parameters:
  simecol::times(model) <- time_seq
  if (!is.null(param)) {
    simecol::parms(model)[names(param)] <- param
  }
  # Define the solver type: see simecol::solver()
  if (!is.null(solver_type)) {
    simecol::solver(model) <- solver_type
  }
  if (is.null(scenarii)) {
    scenarii  <- list(default = simecol::init(model))
  }

  # Define the combination of parameters
  scenar_gradient <- gradient
  scenar_gradient$scenario <- names(scenarii)
  comb <- expand.grid(scenar_gradient) %>%
    dplyr::mutate(
      inits = purrr::map(scenario, function(x) scenarii[[x]])
      )
  param_combination <- dplyr::select(comb, -scenario, -inits) %>%
    df2list()

   run <- parallel::mcMap(
     run_simecol,
     inits = comb[["inits"]],
     param = param_combination,
     MoreArgs = list(model = model)
     )
  # Run the simulations
  #cat(str(run), str(param_combination), str(comb[["inits"]]),
    #sep = "\n")
  output <- as.tibble(comb) %>%
    dplyr::mutate(
    scenario = comb[, "scenario"],
    run = run
    ) %>%
    dplyr::select(-inits) %>%
    dplyr::select(scenario, everything()) %>%
    dplyr::select(-run, everything())

  # Save model parameters
  model <- eval(call(model_spec))
  basis_param <- simecol::parms(model)[which(!names(simecol::parms(model)) %in%
    names(gradient))]
  if (!is.null(param)) {
    basis_param[names(param)] <- param
  }

  return(
    structure(
      list(
	model = model_spec,
	inits = scenarii,
	param = basis_param,
	gradient = gradient,
	run = output),
    class = c("scenarii", "list"))
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

    # Specify the final cover 
    mi_cover <- ini_cover / 2
    mi_low_cover <-  low_cover / 2
    high_cover <- ini_cover - mi_low_cover #For low_P and low_N 

  # three or four states model ?
  variables <- names(simecol::init(model))
  if (all(names(simecol::init(two_facilitation_model())) %in% variables)){
    # four states
    nurse_only <- c(N = ini_cover, P = 0, D = .1,
      NP = 0, PP = 0, NN = ini_cover * ini_cover,
      DD = .1 * .1, PD = 0, ND = ini_cover * .1)
    protegee_only <- c(N = 0, P = ini_cover, D = .1,
      NP = 0, PP = ini_cover * ini_cover, NN = 0,
      DD = .1 * .1, PD = ini_cover * .1, ND = 0)


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

    # The three state model
  } else if (all(names(simecol::init(indirect_facilitation_model())) %in% variables)) {

    nurse_only <- c(N = ini_cover, P = 0,
      NP = 0, PP = 0, NN = ini_cover * ini_cover)
    protegee_only <- c(N = 0, P = ini_cover,
      NP = 0, PP = ini_cover * ini_cover, NN = 0)

    low_N <- c(N = mi_low_cover, P = high_cover,
      NP = high_cover * mi_low_cover, PP = high_cover * high_cover,
      NN = mi_low_cover * mi_low_cover)
    low_P <- c(N = high_cover, P = mi_low_cover,
      NP = high_cover * mi_low_cover, PP = mi_low_cover * mi_low_cover,
      NN = high_cover * high_cover)

    low_together <- c(N = mi_low_cover, P = mi_low_cover,
      NP = mi_low_cover * mi_low_cover, PP = mi_low_cover * mi_low_cover,
      NN = mi_low_cover * mi_low_cover)

    together <- c(N = mi_cover, P = mi_cover,
      NP = mi_cover * mi_cover, PP = mi_cover * mi_cover, NN = mi_cover *
	mi_cover)

  } else {
    stop("The model has not been recognized :/ \n A spelling mistake ?")
  }

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
}

#' Run the model by specifying initial values, parameters and the model  
#' 
#' Run the model over multidimensional gradient and initial values.
#' 
#' @param inits a vector of initial values. See simecole::init 
#' @param params a named vector of length one
#' @param model a function containing a odeModel
#'
#' @return a data.frame 
#' @export
run_simecol <- function(inits, params, model) {

  simecol::parms(model)[names(params)] <- params
  simecol::init(model) <- inits

  run <- simecol::sim(model)
  output <- simecol::out(run) %>%
    dplyr::select(-time)

  return(output)
}
