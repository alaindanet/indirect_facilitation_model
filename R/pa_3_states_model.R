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

    dNP <- PE * Ncolonize(N, NE, E, z, del, b, c, gamma1, nbs = "P") +
      NE * Pcolonize(P, PE, E, z, del, b, c, g, n, nbs = N) -
      2 * NP * die(m)
    dPP <- 2 * PE * Pcolonize(P, PE, E, z, del, b, c, g, n, nbs = "P") -
      2 * PP * die(m)
    dNN <- 2 * NE * Ncolonize(N, NE, E, z, del, b, c, gamma1, nbs = "N") -
      2 * NN * die(m)
    dN <- E * Ncolonize(N, NE, E, z, del, b, c, gamma1, nbs = NULL)  - N * die(m)
    dP <- E * Pcolonize(P, PE, E, z, del, b, c, g, n, nbs = NULL)  - P * die(m)

    list(c(dNP, dPP, dEE, dN, dP))
})
}

#TODO: PB! There is three P colonization equation and three N because of the
#neighbors. Peut - être ajouter des if dans les fonctions ?
