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
    parms = c(z = 4, del = 0.1, b = 0.8, c = 0.2, g = 0.08, m = 0.2, gamma1 =
      0.08, tau = 20),
    init = c(N = .4, P = .4, NP = 0.1, PP = .1, NN = .1),
    solver = "lsoda"
    )
}
