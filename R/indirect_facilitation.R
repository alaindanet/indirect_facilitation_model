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
      s_cost = function(gamma1, k){
	l_cost(gamma1) / (1 + k * l_cost(gamma1)) # Non-linear system
      }
      ),
    times = c(from = 0, to = 100, by = 0.1),
    parms = ,
    init = c(u = 10, v = 5, w = 0.1),
    solver = "lsoda"
    )
}
