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
    } else if (is.character(nbs) & nbs %in% c("P", "N", "D")){
      # Good
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

#' Check integration problem in runs 
#' 
#' Check if there is NaN (Not a Number) or negative densities during the run
#' 
#' @param df a dataframe containing the run of an ODE model
#' @return logical  
#'
#' @export
is_run_normal <- function (df) {
  check <- ifelse(
    any(
      sapply(df, simplify = "matrix", is.nan) |
      df < 0),
    FALSE,
    TRUE)
  return(check)
}

#' Clean averaged simulation runs 
#' 
#' Suppress simulation runs and spurious information
#' @param run a tibble or a data.frame
#' @return a tibble or a data.frame
#' @export
clean_run <- function(x){
  if("PARTITION_ID" %in% names(x)) {
	    x %>% dplyr::select(-runs, -PARTITION_ID)
	  } else {
	    x %>% dplyr::select(-runs)
	  }
}
