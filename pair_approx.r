#'#'#############################################
#'  Pair approximation solutions with deSolve  #'
#'##############################################'

library(deSolve)                       #ODE package
?ode()

FloMod <- function (Time, State, Pars) {  
	with(as.list(c(State, Pars)), {
		     Colonisation <- ( del*rho["+"] + (1-del) * (z-1)/z * (rho["+"] - rho["++"] - rho["+-"])/(1- rho["+"] - rho["-"]) )*(b- c* rho["+"] - g*( 1 - ( z-1 )/z* (rho["+"] - rho["++"] - rho["+-"])/(1- rho["+"] - rho["-"])*n) )
		     Mortality    <- m
		     Degradation  <- d
		     Regeneration <- (r + (z-1)/z*f* (rho["+"] - rho["++"] - rho["+-"])/(1- rho["+"] - rho["-"]))

		     drhopp <- 2 * (rho["+"] - rho["++"] - rho["+-"]) * Colonisation - 2 * rho["++"] * Mortality
		     drhopm <- (rho["+"] - rho["++"] - rho["+-"]) * Degradation + (rho["-"] - rho["--"] - rho["+-"]) * Colonisation - rho["+-"] * (Mortality + Regeneration)
		     drhomm <- 2 * (rho["-"] - rho["--"] - rho["+-"]) * Degradation - 2 * rho["--"] * Regeneration
		     drhom  <- (1- rho["+"] - rho["-"]) * Degradation - rho["-"] * Regeneration
		     drhop <- (1- rho["+"] - rho["-"]) * Colonisation - rho["+"] * Mortality

		     return(list(c(drhopp, drhopm, drhomm, drhom, drhop)))
})
}

pars <- c(del = 0.1,
	  z   = 4,
	  b   = 0.8,
	  c   = 0.2,
	  g   = 0,
	  n   = 1,
	  m   = 0.2,
	  f   = 0.9,
	  d   = 0.1,
	  r   = 0
	  )

# Set the initial values
rhop  <- 0.8
rhom  <- 0.1
rho0  <- 0.1
rho0p <- 0.01
rho0m <- 0.01
rho   <- c("++" = rhop - rho0p - 0.01,
	 "+-"   = 0.01,
	 "--"   = rhom - rho0m - 0.01,
	 "-"    = rhom,
	 "+"    = rhop

	 )

times <- seq(0,500, by=1)

out  <- ode(rho, times, FloMod, pars)
plot(out)

