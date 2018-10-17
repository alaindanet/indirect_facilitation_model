#' Define a custom solver checking for steady state
#' 
#' It stops the numerical integration when the sum of the absolute values of the
#' derivatives become nuls
#' 
#'
#' @export
steady_state_3 <- function(time, init, func, parms) {
  root <- function(time, init, parms) {
    dstate <- unlist(three_states_sys(time, init, parms))
    return(sum(abs(dstate)) - 1e-10)
  }
  lsodar(time, init, func, parms, rootfun = root, atol = 1e-10, rtol = 1e-10)
}
steady_state_4 <- function(time, init, func, parms) {
  root <- function(time, init, parms) {
    dstate <- unlist(four_states_sys(time, init, parms))
    return(sum(abs(dstate)) - 1e-10)
  }
  lsodar(time, init, func, parms, rootfun = root, atol = 1e-10, rtol = 1e-10)
}
