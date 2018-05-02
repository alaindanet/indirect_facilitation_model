#' Colonization function for the nurse pair approximation 
#'
#' 
#' @param N  
#' @param NE  
#' @param E  
#' @param z  
#' @param del  
#' @param b  
#' @param c  
#' @param gamma1  
#' @param nbs  
#'
#' @return a rate of transition
#'
#' @export
Ncolonize <- function(N, NE, E, z, del, b, c, gamma1, nbs){

  NE_proba <- NE_context(nbs, NE, E, z)

  dispersion <- del * N + NE_proba * (1 - del)

  early_survival <- b - c * (1 - E) - gamma1

  return(dispersion * early_survival)

}

#' Colonization function for the nurse pair approximation 
#'
#' 
#' @param P
#' @param PE
#' @param E
#' @param z
#' @param del
#' @param b
#' @param c
#' @param g
#' @param nbs
#'
#' @return a rate of transition
#'
#' @export
Pcolonize <- function(P, N, NE, PE, E, z, del, b, c, g, nbs, p_fun, n, tau){

PE_proba <- PE_context(nbs, PE, E, z)
NE_proba <- NE_context(nbs, NE, E, z)

protection <- p_fun(NE_proba, n, tau)

dispersion <- del * P + PE_proba * (1 - del)

early_survival <- b - c * (1 - E) - g * (1 - protection)

return(dispersion * early_survival)

}

#' Neighboring function 
#'
#' Define the probability of an event according to the Neighboring cells for a 
#' nurse cell. 
#' @param nbs 
#' @param NE 
#' @param E 
#' @param z 
#'
#' @return definition of the probability
#'
#' @export
NE_context <- function(nbs, NE, E, z) {

  check_nbs(nbs)
  check_z(z)

  if (is.null(nbs)){
    return(NE / E)
  } else if (nbs != ("N")) {
    return( ( (z - 1) / z) * NE / E)
  } else {
    return(1 / z + ( (z - 1) / z) * NE / E)
  }

}

#' Neighboring function 
#'
#' Define the probability of an event according to the Neighboring cells for a 
#' protegé cell. 
#' @param nbs 
#' @param NE 
#' @param E 
#' @param z 
#'
#' @return definition of the probability
#'
#' @export
PE_context <- function(nbs, PE, E, z) {

  check_nbs(nbs)
  check_z(z)

  if (is.null(nbs)){
    return(PE / E)
  } else if (nbs != "P") {
    return( ( (z - 1) / z) * PE / E)
  } else {
    return(1 / z + ( (z - 1) / z) * PE / E)
  }

}

#' Neighboring function 
#'
#' Define the probability of an event according to the Neighboring cells for a 
#' degraded cell. 
#' @param nbs 
#' @param ND
#' @param PD
#' @param D 
#' @param z 
#'
#' @return definition of the probability
#'
#' @export
D_context <- function(nbs, ND, PD, D, z) {

  check_nbs(nbs)
  check_z(z)

  if (is.null(nbs)){
    return(ND / D)
  } else if (nbs != "D") {
    return( 1 / z + ( (z - 1) / z) * (ND / D + PD / D))
  } else {
    return( ( (z - 1) / z ) * (ND / D + PD / D))
  }

}

#' Death transition function 
#'
#' Define the probability to die for a cell.
#'
#' @param m rate of mortality.  1/m defines the life expectancy.  
#'
#' @return definition of the probability
#'
#' @export
die <- function(m) {
  return(m)
}

#' Degradation transition function 
#'
#' Define the probability to be degraded for a cell.
#'
#' @param d rate of degradation.  1/d defines the degradation time expectancy.  
#'
#' @return definition of the probability
#'
#' @export
degradation <- function(d) {
  return(d)
}

#' Degradation transition function 
#'
#' Define the probability to be degraded for a empty cell.
#'
#' @param d rate of degradation.  1/d defines the degradation time expectancy.  
#'
#' @return definition of the probability
#'
#' @export
degrade <- function(d) {
  return(d)
}

#' Regeneration transition function 
#'
#' Define the probability to be regenerated for a degraded  cell.
#'
#' @param r rate of regeneration.  1/r defines the expected time before
#' regeneration.
#'
#' @return definition of the probability
#'
#' @export
regen <- function(r) {
  return(r)
}

#' Facilitation transition function 
#'
#' Define the probability to be regenerated for a degraded cell thanks to the
#' improvement of local conditions by an adult plant.
#'
#' @param f rate of local condition improvement.
#'
#' @return definition of the probability
#'
#' @export
facilitate <- function(ND, PD, D, z, f, nbs = NULL) {

  D_proba <- D_context(nbs, ND, PD, D, z)
  facilitation <- f * D_proba

  return(facilitation)
}

#' Protection function for the protegee 
#' 
#' @param type  
#' @param NE_proba 
#'
#' @return a rate of transition
#'
#' @export
compute_as <- function(type = "linear"){

  if (is.null(type)) {
    # No protection
    protection <- function (NE_proba, n, tau) {
      return(0)
    }

  } else if (type == "linear") {
    protection <- function (NE_proba, n, tau) {
      return(NE_proba * n)
    }

  } else if (type == "first_protect") {
    protection <- function (NE_proba, n, tau) {
      return(1 - exp(- tau * NE_proba))
    }
  }
  return(protection)
}

#' relationship between the effect of a nurse and the cost of defense 
#' 
#' @param gamma1 cost of defence
#' @param u shape of the curve.  
#'
#' @return a rate of transition
#'
#' @export
compute_p_one_z <- function(gamma1, u) {
  return(1 - exp(- u * gamma1))
}

#' Define the shape of the relationship between the protection effect and the
#' number of nurse 
#' 
#' @param p the effect of one nurse 
#' @param z the number of neighbors
#'
#' @return tau 
#'
#' @export
compute_tau <- function(p = 9 / 10, z = 4) {
  log(-1 / (p - 1)) / (1 / z)
}
