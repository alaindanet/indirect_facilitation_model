#' pair approximation system of ODEs with a nurse, a protege and empty states
#'
#' 
#' @param ... Nothing 
#' @return A list of ODEs 
#'
#' @export

three_states_sys <- function(time, init, parms) {
  # Variables:
  NP <- init["NP"]
  NN <- init["NN"]
  PP <- init["PP"]
  N <- init["N"]
  P <- init["P"]
  # Variables coming for conservation law:
  PE <- P - NP - PP
  NE <- N - NP - NN
  E <- 1 - N - P

  with(as.list(parms), {
    # Compute the associative protection
    n <- asymp_cost(gamma1, tau)
    #Thresholds
    #N <- ifelse(N < extinction_threshold, 0, N)
    #P <- ifelse(P < extinction_threshold, 0, P)

    dNP <- PE * Ncolonize(N, NE, E, z, del, b, c, gamma1, nbs = "P") +
      NE * Pcolonize(P, N, NE, PE, E, z, del, b, c, g, n, nbs = "N") -
      2 * NP * die(m)
    dPP <- 2 * PE * Pcolonize(P, N, NE, PE, E, z, del, b, c, g, n, nbs = "P")-
      2 * PP * die(m)
    dNN <- 2 * NE * Ncolonize(N, NE, E, z, del, b, c, gamma1, nbs = "N") -
      2 * NN * die(m)
    dN <- E * Ncolonize(N, NE, E, z, del, b, c, gamma1, nbs = NULL) -
      N * die(m)
    dP <- E * Pcolonize(P, N, NE, PE, E, z, del, b, c, g, n, nbs = NULL) -
      P * die(m)

    # the variables should be returned in the same order the init values
    # init = c(N = .4, P = .4, NP = 0.3, PP = .3, NN = .3),
    list(c(dN, dP, dNP, dPP, dNN))
})
}

