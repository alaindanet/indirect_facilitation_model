#' Define a custom solver checking for steady state
#' 
#' It stops the numerical integration when the sum of the absolute values of the
#' derivatives become nuls
#' 
#'
#' @export
steady_state_3 <- function(time, init, func, parms) {
  root <- function(time, init, parms) {
    dstate <- unlist(three_states_sys(time, init, parms))
    return(sum(abs(dstate)) - 1e-10)
  }
  lsodar(time, init, func, parms, rootfun = root, atol = 1e-10, rtol = 1e-10)
}
steady_state_4 <- function(time, init, func, parms) {
  root <- function(time, init, parms) {
    dstate <- unlist(four_states_sys(time, init, parms))
    return(sum(abs(dstate)) - 1e-10)
  }
  lsodar(time, init, func, parms, rootfun = root, atol = 1e-10, rtol = 1e-10)
}

ca_solver <- function(y, times=NULL, func=NULL, parms=NULL,
  animate=FALSE, ...) {
  observer <- function(landscape) {
    # nurse, protegee, empty, degraded
    rho_nurse    <- sum(landscape == 1) / length(landscape)
    rho_protegee <- sum(landscape == 2) / length(landscape)
    rho_empty    <- sum(landscape == 3) / length(landscape)
    rho_degraded <- sum(landscape == 4) / length(landscape)
    neigh_n      <- fourneighbors(landscape, state = 1, bounds = 1)
    neigh_p      <- fourneighbors(landscape, state = 2, bounds = 1)
    neigh_veg    <- neigh_n + neigh_p

    qnp  <- mean(neigh_n[landscape == 2])
    qpn  <- mean(neigh_p[landscape == 1])
    qnn  <- mean(neigh_n[landscape == 1])
    qpp  <- mean(neigh_p[landscape == 2])
    qveg <- mean(neigh_veg[landscape %in% c(1, 2)])

    c(
      N    = rho_nurse,
      P = rho_protegee,
      E    = rho_empty,
      D = rho_degraded,
      qnp      = qnp,
      qpn      = qpn,
      qnn      = qnn,
      qpp      = qpp,
      qveg     = qveg
      )
  }
  init <- y@init
  times <- fromtoby(y@times)
  func <- y@main
  parms <- y@parms
  inputs <- y@inputs
  equations <- y@equations
  equations <- addtoenv(equations)
  environment(func) <- environment()
  parms$DELTAT <- 0
  # Check config
  unstable <- TRUE
  two_species <- TRUE
  nb_check <- 0
  check_points <- seq(length(times) %/% 200) * 200
  i <- 1
  out <- as.list(seq_along(times))
  res <- observer(init)
  out[[i]] <- res
  # Average densities 
  no_val <- rep(NA, length(check_points))
  avg <- tibble::tibble(N = no_val, P = no_val)

  # Loop:
  while (i < length(times) & unstable & two_species) {
    i <- i + 1
    time <- times[i]
    parms$DELTAT <- times[i] - times[i-1]
    init <- func(time, init, parms)
    res <- observer(init)
    out[[i]] <- res
    
    # Conditions check  
    if (i >= 200 & i %in% check_points) {

      nb_check <- nb_check + 1
      # Select data:
      test <- as.data.frame(do.call(rbind, out[1:i]))
      if (nb_check == 1){
	avg[nb_check, c("N", "P")] <- sapply(test[1:check_points[nb_check], c("N", "P")], mean)
      } else {
	avg[nb_check, c("N", "P")] <- sapply(test[
	  check_points[nb_check - 1]:check_points[nb_check],
	  c("N", "P")], mean)
	# Test for stability
	stab_stat <- sapply(avg, function(x) {
	  abs(x[nb_check - 1] - x[nb_check]) / x[nb_check - 1]
	  })
	if (all(stab_stat < .005)) {
	  unstable <- FALSE
	}
      }
      # Test for presence of the two species:
      two_species_stat <- sapply(avg, function(x) any(x[nb_check] == 0))
      if (any(two_species_stat)){
	two_species <- FALSE
      }
    }
  }
  out <- do.call(rbind, out[1:i])
  row.names(out) <- NULL
  out <- cbind(time = times[1:i], out)
  as.data.frame(out)
}
