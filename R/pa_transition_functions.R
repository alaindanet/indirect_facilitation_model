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
Pcolonize <- function(P, PE, E, z, del, b, c, g, nbs){

  PE_proba <- PE_context(nbs, PE, E, z)

  dispersion <- del * P + PE_proba * (1 - del)

  early_survival <- b - c * (1 - E) - g

  return(dispersion * early_survival)

}

#' Neighboring function 
#'
#' 
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
    return(NE)
  } else if (nbs == "P") {
    return( ( (z - 1) / z) * NE / E)
  } else {
    return(1 / z + ( (z - 1) / z) * NE / E)
  }

}

#' Neighboring function 
#'
#' 
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
    return(PE)
  } else if (nbs == "N") {
    return( ( (z - 1) / z) * PE / E)
  } else {
    return(1 / z + ( (z - 1) / z) * PE / E)
  }

}
