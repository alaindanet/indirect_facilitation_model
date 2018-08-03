four_states_ca <- function(time, init, parms) {

  # Variables:
  landscape <- init
  rho_nurse <- sum(landscape == 1) / length(landscape)
  rho_protegee <- sum(landscape == 2) / length(landscape)

  # Neighbors:
  if(parms$z != 4) {
    stop("z should be equal to 4")
  }
  neigh_n <- fourneighbors(landscape, state = 1, bounds = 1)
  neigh_p <- fourneighbors(landscape, state = 2, bounds = 1)

  # Compute the effect of a nurse (used if protection type if protection_type is
  # first_protect)
    if (unlist(parms["protection_type"]) == "first_protect") {
      p_one_z <- compute_p_one_z(unlist(parms["gamma1"]), unlist(parms["u"]))
      parms$tau <- compute_tau(p_one_z, unlist(parms["z"]))
    } else {
      stop("only the 'first protect' protection type is implemented")
    }


  with(parms, {

    colonization_nurse <- (del * rho_nurse + (1 - del) * neigh_n / z) *
      (b - c * (rho_protegee + rho_nurse) - gamma1) * DELTAT

    colonization_protegee <- (del * rho_protegee + (1 - del) * neigh_p / z) *
      (b - c * (rho_protegee + rho_nurse) - g * (1 - (1 - exp(- tau * neigh_n / z)))) * DELTAT

    #cat(
      #"protegee: ", mean(colonization_protegee), "\n",
      #"nurse: ", mean(colonization_nurse), "\n",
      #"grazing effect: ", mean(g * (1 - (1 - exp(- tau * neigh_n / z)))), "\n", #g *
      #"neigh_n: ", mean(neigh_n / z), "\n",
      #"neigh_p: ", mean(neigh_p / z ), "\n"
      #)
    death <- m * DELTAT
    # calculate regeneration rate and degradation rate
    regeneration <- (r + f * (neigh_n + neigh_p) / z) * DELTAT
    degradation <- d * DELTAT
    # Apply rules
    rnum <- runif(length(landscape)) # one random number between 0 and 1 for each cell
    new_landscape <- landscape

    ## New nurses
    new_landscape[which(landscape == 3 & rnum <= colonization_nurse)] <- 1

    ## New protegees
    new_landscape[which(landscape == 3 & rnum > colonization_nurse &
      rnum <= colonization_protegee + colonization_nurse)] <- 2

    ## New degraded cells
    new_landscape[which(landscape == 3 & rnum > colonization_nurse +
      colonization_protegee & rnum <= colonization_nurse +
      colonization_protegee + degradation)] <- 4

    ## New empty
    new_landscape[which(landscape == 1 & rnum <= death)] <- 3
    new_landscape[which(landscape == 2 & rnum <= death)] <- 3
    new_landscape[which(landscape == 4 & rnum <= regeneration)] <- 3

    # check for sum of probabilities to be inferior 1 and superior 0
    if (any(c(colonization_nurse + degradation,
	  colonization_protegee + degradation,
	  colonization_nurse + colonization_protegee + degradation,
	  death, regeneration) > 1 )) {
      warning("a set probability is exceeding 1 in run! decrease delta!!!")
    }
    if (any(c(colonization_nurse, colonization_protegee,
	  degradation, regeneration,
	  death) < 0)) {
      warning("a set probability falls below 0 in run balance parameters!!!")
    }

    new_landscape

})
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
  check_points <- seq(length(times) %/% 300) * 300
  i <- 1
  out <- as.list(seq_along(times))
  res <- observer(init)
  out[[i]] <- res
  while (i < length(times) & unstable & two_species) {
    # Loop:
    i <- i + 1
    time <- times[i]
    parms$DELTAT <- times[i] - times[i-1]
    init <- func(time, init, parms)
    res <- observer(init)
    out[[i]] <- res

    # Conditions check  
    if (i >= 300 & i %in% check_points) {

      nb_check <- nb_check + 1
      # Select data:
      test <- do.call(rbind, out[1:i])
      if(nb_check == 1){
	checked_data <- test[1:check_points[nb_check], c("N", "P")]
      } else {
	checked_data <- test[
	  check_points[nb_check - 1]:check_points[nb_check],
	  c("N", "P")]
      }
      # Measurements:
      test <- sapply(as.data.frame(checked_data), function(x) all(x == 0))
      if (any(test)){
	two_species <- FALSE
      }
      if (two_species) { # If one species is extinct, the skewness returns NaN
	test <- sapply(as.data.frame(checked_data), moments::skewness)
	if (all(test < .25)) {
	  unstable <- FALSE
	}
      }
    }
  }
  out <- do.call(rbind, out[1:i])
  row.names(out) <- NULL
  out <- cbind(time = times[1:i], out)
  as.data.frame(out)
}
