#' pair approximation system of ODEs with a nurse, a protege, empty and degraded states
#'
#' 
#' @param ... Nothing 
#' @return A list of ODEs 
#'
#' @export
four_states_sys <- function(time, init, parms) {
  # Variables:
  NN <- init["NN"]
  PP <- init["PP"]
  DD <- init["DD"]

  ND <- init["ND"]
  PD <- init["PD"]
  NP <- init["NP"]

  N <- init["N"]
  P <- init["P"]
  D <- init["D"]

  # Variables coming for conservation law:
  PE <- P - NP - PP - PD
  NE <- N - NP - NN - ND
  DE <- D - DD - PD - ND
  E <- 1 - N - P - D

  # Define the associative protection function
  p_fun <- compute_as(type = unlist(parms["protection_type"]))

  # Compute the effect of a nurse (used if protection type if protection_type is
  # first_protect)
  if (!is.null(unlist(parms["protection_type"]))) {
    if (unlist(parms["protection_type"]) == "first_protect") {
      p_one_z <- compute_p_one_z(unlist(parms["gamma1"]), unlist(parms["u"]))
      tau <- compute_tau(p_one_z, unlist(parms["z"]))
    } else if (unlist(parms["protection_type"]) == "linear") {
      tau <- unlist(parms["tau_n"])
      parms["n"] <- list(1 - exp(- tau * unlist(parms["gamma1"])))
    }
  }

  with(as.list(parms), {
    #Thresholds
    #N <- ifelse(N < extinction_threshold, 0, N)
    #P <- ifelse(P < extinction_threshold, 0, P)

    dNP <- PE * Ncolonize(N, NE, E, z, del, b, c, gamma1, nbs = "P") +
      NE * Pcolonize(P, N, NE, PE, E, z, del, b, c, g, nbs = "N", p_fun, n, tau) -
      2 * NP * die(m)
    dND <- NE * degrade(d) +
      DE * Ncolonize(N, NE, E, z, del, b, c, gamma1, nbs = "D") -
      ND * die(m) -
      ND * (regen(r) + facilitate(ND, PD, D, z, f, nbs = "N"))
    dPD <- PE * degrade(d) +
      DE * Pcolonize(P, N, NE, PE, E, z, del, b, c, g, nbs = "D", p_fun, n, tau) -
      PD * die(m) -
      PD * (regen(r) + facilitate(ND, PD, D, z, f, nbs = "P"))

    dPP <- 2 * PE * Pcolonize(P, N, NE, PE, E, z, del, b, c, g, nbs = "P", p_fun, n, tau) -
      2 * PP * die(m)
    dNN <- 2 * NE * Ncolonize(N, NE, E, z, del, b, c, gamma1, nbs = "N") -
      2 * NN * die(m)
    dDD <- 2 * DE * degrade(d) -
      2 * DD * (regen(r) + facilitate(ND, PD, D, z, f, nbs = NULL))

    # Singleton variables
    dN <- E * Ncolonize(N, NE, E, z, del, b, c, gamma1, nbs = NULL) -
      N * die(m)
    dP <- E * Pcolonize(P, N, NE, PE, E, z, del, b, c, g, nbs = NULL, p_fun, n, tau) -
      P * die(m)
    dD <- E * degrade(d) - D * (regen(r) + facilitate(ND, PD, D, z, f, nbs = NULL))

    # the variables should be returned in the same order the init values
    # init = c(N = .4, P = .4, NP = 0.3, PP = .3, NN = .3),
    list(c(dN, dP, dD, dNP, dPP, dNN, dDD, dPD, dND))
})
}
