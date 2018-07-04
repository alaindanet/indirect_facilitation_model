
#' Convert a data.frame to a named list
#'
#' @param x a data.frame 
#' @return a list. 
#' @details Each element of the list contains a named vector. The names
#' are those of the columns of the input.
#'
#' @export
df2list <- function (x) {
  test <- apply(x, 1, as.list)
  lapply(test, unlist)
}
