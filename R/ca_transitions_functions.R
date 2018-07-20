#' Compute number of neighbors out of 4 neighboring cells 
#'
#' This function is based on the simecol::neighbors function.
#' @inherit simecol::neighbors 
#'
#' @return a matrix
#' @export
fourneighbors <- function(landscape, state = 1, bounds = 1) {

  neighborhood <- matrix(
    c(0, 1, 0,
      1, 0, 1,
      0, 1, 0),
    nrow = 3)
  nb_neighbors <- simecol::neighbors(
    x = landscape,
    state = state,
    wdist = neighborhood,
    bounds = bounds)

}
