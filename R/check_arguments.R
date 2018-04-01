#' Check the neighbor arguments 
#'
#' \code{indirect_facilitation_model} returns the simObj.
#' 
#' @param nbs string. Either "N", "P" or "NULL"
#' @return stop if the argument is badly define 
#'
#' @export
check_nbs <- function(nbs) {
  if (is.null(nbs)) {
      # Good
    } else if (is.character(nbs) & nbs %in% c("P", "N")){
    } else {
    stop("nbs is badly defined")
  }
}

#' Check the neighbor arguments 
#' 
#' 
#' @param z numeric. Either 4 or 8 
#' @return stop if the argument is badly defined 
#'
#' @export
check_z <- function(z) {
  if (is.numeric(z)) {

    if (z %in% c(4, 8)){
      # Good
    } else {
      stop("z is badly defined")
    }
  } else {
    stop("z is badly defined")
  }
}

#' Check the neighbor arguments 
#' 
#' 
#' @param arg vector List of arguments 
#' @return stop if the argument is badly defined 
#'
#' @export
check_var <- function(arg) {

  !is.vector(arg) || stop("arg should be a vector")
  !is.numeric(arg) || stop("variables should be numeric")

}
