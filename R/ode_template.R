#' Template to bluid an ode system
#'
#' \code{upca} returns the .
#' 
#' This a template from the simecol package. It produces guidelines to easily
#' compute odes.
#'
#' @param ... Numeric, complex, or logical vectors.
#' @param na.rm A logical scalar. Should missing values (including NaN)
#'   be removed?
#' @return If all inputs are integer and logical, then the output
#'   will be an integer. If integer overflow
#'   \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
#'   will be NA with a warning. Otherwise it will be a length-one numeric or
#'   complex vector.
#'
#'   This example is coming from the documentation of the simecol package. See
#'   \url{http://simecol.r-forge.r-project.org/} for more details.
#' @examples
#' sum(1:10)
#' sum(1:5, 6:10)
#' sum(F, F, F, T, T)
#'
#' sum(.Machine$integer.max, 1L)
#' sum(.Machine$integer.max, 1)
#'
#' \dontrun{
#' sum("a")
#' }

upca <- new("odeModel",
  main = function(time, init, parms) {
    u <- init[1]
    v <- init[2]
    w <- init[3]
    with(as.list(parms), {
      du <- a * u
      - alpha1 * f(u, v, k1)
      dv <- -b * v
      + alpha1 * f(u, v, k1) - alpha2 * f(v, w, k2)
      dw <- -c * (w - wstar) + alpha2 * f(v, w, k2)
      list(c(du, dv, dw))
})
  },
  equations = list(
    f1 = function(x, y, k){x * y}, # Lotka-Volterra
    f2 = function(x, y, k){f1(x, y, k) / (1 + k * x)} # Holling II
    ),
  times = c(from = 0, to = 100, by = 0.1),
  parms = c(a = 1, b = 1, c = 10, alpha1 = 0.2, alpha2 = 1,
    k1 = 0.05, k2 = 0, wstar = 0.006),
  init = c(u = 10, v = 5, w = 0.1),
  solver = "lsoda"
  )
equations(upca)$f <- equations(upca)$f1
