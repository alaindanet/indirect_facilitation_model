
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
  time_seq = NULL,
  param = NULL, nb_cores = NULL, solver_type = NULL,
  scenarii = NULL, set_tail = NULL, nrep = NULL) {

  if (is.null(gradient)) {
    stop("Please provide a gradient")
  }

  model <- eval(call(model_spec))
  # Define parameters:

  if (!is.null(time_seq)) {
  simecol::times(model) <- time_seq
  }
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
  if (!is.null(nrep)) {
  gradient$rep <- seq.int(1, nrep)
  }
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
     MoreArgs = list(model = model, set_tail = set_tail)
     )
  # Run the simulations
  output <- as.tibble(comb) %>%
    dplyr::mutate(
    scenario = comb[, "scenario"],
    run = run
    ) %>%
    dplyr::select(-inits) %>%
    dplyr::select(scenario, dplyr::everything()) %>%
    dplyr::select(-run, dplyr::everything())

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
	model = eval(call(model_spec)),
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
#' @param type character. the scenario. Can be either "all", "bifurcation",
#' "together", "low_together", "nurse_only", "protegee_only", "low_N", "low_P"
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
  ini_cover = .8, low_cover = .001) {

  stopifnot(type %in% c("nurse", "protegee", "together", "low_N", "low_P",
      "low_together", "low_nurse_only", "low_protegee_only", "all", "bifurcation",
      "nurse_bifurcation", "protegee_bifurcation"))

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
    low_nurse_only <- c(N = low_cover, P = 0, D = .99,
      NP = 0, PP = 0, NN = low_cover * low_cover,
      DD = .99 * .99, PD = 0, ND = low_cover * .99)
    protegee_only <- c(N = 0, P = ini_cover, D = .1,
      NP = 0, PP = ini_cover * ini_cover, NN = 0,
      DD = .1 * .1, PD = ini_cover * .1, ND = 0)
    low_protegee_only <- c(N = 0, P = low_cover, D = .99,
      NP = 0, PP = low_cover * low_cover, NN = 0,
      DD = .99 * .99, PD = low_cover * .99, ND = 0)

    low_N <- c(N = mi_low_cover, P = high_cover, D = .1,
      NP = high_cover * mi_low_cover, PP = high_cover * high_cover,
      NN = mi_low_cover * mi_low_cover, DD = .1 * .1,
      PD = high_cover * .1, ND = mi_low_cover * .1)
    low_P <- c(N = high_cover, P = mi_low_cover, D = .1,
      NP = high_cover * mi_low_cover, PP = mi_low_cover * mi_low_cover,
      NN = high_cover * high_cover, DD = .1 * .1,
      PD = mi_low_cover * .1, ND = high_cover * .1)

    low_together <- c(N = mi_low_cover, P = mi_low_cover, D = .99,
      NP = mi_low_cover * mi_low_cover, PP = mi_low_cover * mi_low_cover,
      NN = mi_low_cover * mi_low_cover, DD = .99 * .99,
      PD = mi_low_cover * .99, ND = mi_low_cover * .99)

    together <- c(N = mi_cover, P = mi_cover, D = .1,
      NP = mi_cover * mi_cover, PP = mi_cover * mi_cover,
      NN = mi_cover * mi_cover, DD = .1 * .1,
      PD = mi_cover * .1, ND = mi_cover * .1)

    # The three state model
  } else if (all(names(simecol::init(indirect_facilitation_model())) %in% variables)) {

    nurse_only <- c(N = ini_cover, P = 0,
      NP = 0, PP = 0, NN = ini_cover * ini_cover)
    low_nurse_only <- c(N = low_cover, P = 0,
      NP = 0, PP = 0, NN = low_cover* low_cover)
    protegee_only <- c(N = 0, P = ini_cover,
      NP = 0, PP = ini_cover * ini_cover, NN = 0)
    low_protegee_only <- c(N = 0, P = low_cover,
      NP = 0, PP =  low_cover* low_cover, NN = 0)

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
      "low_protegee_only" = low_protegee_only,
      "nurse_only" = nurse_only,
      "low_nurse_only" = low_nurse_only,
      "together" = together,
      "low_N" = low_N,
      "low_P" = low_P,
      "low_together" = low_together
      )
    if(all(type == "all")) {
      return(all_inits)
    } else if (all(type %in% c("protegee_bifurcation", "bifurcation"))) {
      return(all_inits[c("low_protegee_only", "protegee_only", "low_together", "together")])
    } else if (all(type == "bifurcation")) {
      return(all_inits[c("low_together", "together")])
    } else if (all(type == "protegee_bifurcation")) {
      return(all_inits[c("low_protegee_only", "protegee_only")])
    } else if (all(type == "nurse_bifurcation")) {
      return(all_inits[c("low_nurse_only", "nurse_only")])
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
run_simecol <- function(inits, params, model, set_tail = NULL) {

  simecol::parms(model)[names(params)] <- params
  simecol::init(model) <- inits

  run <- simecol::sim(model)

  # For big simulations, keep only the last step
  if (!is.null(set_tail)) {
    output <- simecol::out(run) %>%
      tail(., set_tail)
  } else {
    output <- simecol::out(run)
  }
  #output %<>% select(-time)

  return(output)
}
