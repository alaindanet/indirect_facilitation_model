#' Indirect facilitation model with a nurse and a protege 
#'
#' \code{indirect_facilitation_model} returns the simObj.
#' 
#' @param ... Nothing 
#' @return An ecological model of the class simObj with default definition.
#'
#' @examples
#' model <- indirect_facilitation_model() 
#'
#' \dontrun{
#'  out <- sim(model)
#' 
#' }
#'
#' @export
indirect_facilitation_model <- function() {
  new("odeModel",
    main = three_states_sys,
    equations = list(
      l_cost = function(gamma1){ # linear cost
	gamma1
      },
      asymp_cost = function(gamma1, tau){ return(1 - exp(tau * gamma1))}
      ),
    times = c(from = 0, to = 100, by = 0.1),
    parms = c(z = 4, del = 0.1, b = 0.8, c = 0.2, g = 0.08, m = 0.2, gamma1 = 0.08,
      extinction_threshold = 1 * 10 ^ -3, protection_type = list("linear"), u = 15, n = 1, tau_n = 15),
    init = c(N = .4, P = .4, NP = 0.1, PP = .1, NN = .1),
    solver = "lsoda"
    )
}

two_facilitation_model <- function() {
  new("odeModel",
    main = four_states_sys,
    times = c(from = 0, to = 100, by = 0.1),
    parms = c(z = 4, del = .1, b = .8, c = .2, g = .08, m = .2,
      gamma1 = .08, r = .01, f = .9, d = .1,
      extinction_threshold = 1 * 10 ^ -3, protection_type = list("first_protect"),
      u = 15, n = 1, tau_n = 15),
    init = c(N = .4, P = .4, D = .1,
      NP = .1, PP = .1, NN = .1,
      DD = .1, PD = .1, ND = .1),
    solver = "lsoda"
    )
}

ca_two_facilitation_model <- function () {
  new("gridModel",
    main = four_states_ca,
    times = c(from = 0, to = 100, by = 1),
    parms = c(z = 4, del = .1, b = .8, c = .2, g = .1, m = .2,
      gamma1 = .1, r = .01, f = .9, d = .1, protection_type = list("first_protect"),
      u = 0, skew_threshold = .5),
    init = matrix(sample.int(4, size = 100*100, replace = TRUE, prob = c(.4, .4, .1, .1)), nrow = 100, ncol = 100),
    solver = "ca_solver"
    
    )
}
