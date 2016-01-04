
#'#'####################################################
#'  Four states equations - pair approximation model  #'
#'#####################################################'

library(deSolve)

# TODO: Add the new variables.
# TODO: Define the replacement rules.
# TODO: Write the new equations.
# TODO: define a rule to signal the cases in which the probabilities fall below 0.

#'#'################################################################################################################
#'  What works :
# - Replace the numeric indentation of rho by the character indentation.
# - Defining the replacement rules insides the differential equations and call them after. 

#' What don't:
# - Define the replacement rules in a list and call it in each differential equation functions.
#'################################'


rho <- c(rhonp = rhoini*0.4, # rho[1]
	 rhopp = rhoini*0.4 , # rho[2]
	 rhonn = (1-rhoini)*0.6, # rho[3]
	 rhop  = rhoini,  # rho[4]
	 rhon  = (1-rhoini)*0.9   # rho[5]
		          )# }}}
#
replacement <- list( rhono = rho[5] - rho[3] - rho[1] - rho[8],
	rhopo = ( rho[4] - rho[2] - rho[1] - rho[9] ),
	rhoom = ( rho[6] - rho[7] - rho[8] - rho[9] ),
	rhoo = ( 1 - rho[4] - rho[5] - rho[6] ),
	Qno = rhono/rhoo,
	Qnp = rho[1]/ rho[4],
	Qpn = rho[1]/ rho[5]
)


# Differential equations# {{{

drhonp <- function(rho, parms, z = 4, replacement) { 
rhopo = (rho[4]- rho[2] - rho[1])

res <-   with(parms,
rhopo * (del * rho["rhon"] + (1-del) * (z-1)/z * (rho[5] - rho["rhonn"] - rho[1])/(1- rho[5] - rho[4]))*(b- cn* (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4]) - cpn* (rho[4] - rho[2] - rho[1])/(1- rho[5] - rho[4]) ) + (rho[5] - rho[3] - rho[1]) * (del * rho[4] + (1-del) * (z-1)/z * (rho[4] - rho[2] - rho[1])/(1- rho[5] - rho[4])) * (b- cp* (rho[4] - rho[2] - rho[1])/(1- rho[5] - rho[4]) - cnp* (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4]) - g*( 1 - 1/z*n - ( z-1 )/z * (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4])*n) ) - 2 * rho[1] * m

		          )
return(res)
}

drhopp <- function(rho, parms, z = 4) { 
	  with(parms,
2 * (rho[4] - rho[2] - rho[1]) *
		     (del * rho[4] + (1-del)/z +  (1-del) * (z-1)/z * (rho[4] - rho[2] - rho[1])/(1- rho[5] - rho[4])) *
		     (b- cp* (rho[4] - rho[2] - rho[1])/(1- rho[5] - rho[4]) - cnp* (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4]) - g*( 1 - ( z-1 )/z * (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4])*n) ) -
		     2 * rho[2] * m


	          )
}

drhonn <- function(rho, parms, z = 4) { 
	  with(parms,
2 * (rho[5] - rho[3] - rho[1]) *
		     ( del*rho[5] + (1-del)/z + (1-del) * (z-1)/z * (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4]) )*
		     (b- cn* (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4]) - cpn* (rho[4] - rho[2] - rho[1])/(1- rho[5] - rho[4]) ) - 2 * rho[3] * m  
	       )
}

drhop <- function(rho, parms, z = 4) { 
	  with(parms,                             
	(1- rho[5] - rho[4]) * (del * rho[4] +  (1-del) * (rho[4] - rho[2] - rho[1])/(1- rho[5] - rho[4])) *
		     (b- cp* (rho[4] - rho[2] - rho[1])/(1- rho[5] - rho[4]) - cnp* (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4]) - g*( 1 - (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4])*n) ) - rho[4] * m 
		             )
}

drhon <- function(rho, parms, z = 4) { 
	  with(parms,
(1- rho[4] - rho[5]) * ( del*rho[5] + (1-del) * (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4]) ) *
		     (b- cn* (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4]) - cpn* (rho[4] - rho[2] - rho[1])/(1- rho[5] - rho[4]) ) - rho[5] * m  
		           )
}

odesys <- function (t, rho, parms = modelparms) {
	  list(c(
		     drhonp(rho , parms) ,
		     drhopp(rho , parms) ,
		     drhonn(rho , parms) ,
		     drhop(rho  , parms) ,
		     drhon(rho  , parms)
		         ))
}
# }}}
# Parameters definition# {{{
modelparms <- list(
		   del = 0.1,
		   b   = 0.8,
		   cn  = 0.2,
		   cp  = 0.2,
		   cnp = 0.05,
		   cpn = 0.2,
		   g   = 0.19,
		   n   = 0,
		   m   = 0.2
		   )

rhoini <- 0.1

rho <- c(rhonp = rhoini*0.4, # rho[1]
	 rhopp = rhoini*0.4 , # rho[2]
	 rhonn = (1-rhoini)*0.6, # rho[3]
	 rhop  = rhoini,  # rho[4]
	 rhon  = (1-rhoini)*0.9   # rho[5]
		          )# }}}
# 	 rhom  = 0,                  # rho[6]
# 	 rhomm = 0,                  # rho[7]
# 	 rhonm = 0,                  # rho[8] 
# 	 rhopm = 0                   # rho[9] 

# Replacement rules:

#rhono = rho[5] - rho[3] - rho[1] - rho[8]
#rhopo = rho[4] - rho[2] - rho[1] - rho[9]
#rhoom = rho[6] - rho[7] - rho[8] - rho[9]
#rhoo = 1 - rho[4] - rho[5] - rho[6]
#Qno = rhono/rhoo
#Qnp = rho[1]/ rho[4]
#Qpn = rho[1]/ rho[5]

# }}}
# Model run and plot output# {{{
runmodel <- ode(y = rho, func = odesys, times = seq(1, 300, 30), parms = modelparms)

# transfer into ouput and calculate missing rho values
out       <- as.data.frame(runmodel)
out  <- round(out,3)
out$rho0  <- 1 - out$rhop - out$rhon
out$rhop0 <- with(out, rhop-rhonp-rhopp)
out$rhon0 <- with(out, rhon-rhonp-rhonn)
out$rho00 <- with(out, rho0-rhon0-rhop0)
out$Cnp   <- with(out, (rhonp/rhon)/rhop)

# out$Cpn <- with(out, rhonp
# plot
plot(Cnp ~ time, data = out, type = "l", lwd = 2, ylab = "Species clustering", ylim = c(0,2))
lines(rhop ~ time, data = out, lty = 3, col="red")
lines(rhon ~ time, data = out, lty = 2, col="green")
# }}}
