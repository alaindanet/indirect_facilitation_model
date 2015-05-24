#'#'##################################################################################
#'  Simulation pair-approximation model for two species and associative resistance
#' Author: A Danet
#'###################################################################################'

rm(list=ls())

#
setwd("/home/alain/")


# Equations différentielles du système# {{{
drhonp <- function(rho, parms, z = 4) { 
	  with(parms,
(rho[4]- rho[2] - rho[1]) * (del * rho[5] + (1-del) * (z-1)/z * (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4]))*(b- cn* rho[5] - cpn* rho[4] ) + (rho[5] - rho[3] - rho[1]) * (del * rho[4] + (1-del) * (z-1)/z * (rho[4] - rho[2] - rho[1])/(1- rho[5] - rho[4])) * (b- cp* rho[4] - cnp* rho[5] - g*( 1 - 1/z*n - ( z-1 )/z * (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4])*n) ) - 2 * rho[1] * m

		          )
}

drhopp <- function(rho, parms, z = 4) { 
	  with(parms,
2 * (rho[4] - rho[2] - rho[1]) *
		     (del * rho[4] + (1-del)/z +  (1-del) * (z-1)/z * (rho[4] - rho[2] - rho[1])/(1- rho[5] - rho[4])) *
		     (b- cp* rho[4] - cnp* rho[5] - g*( 1 - ( z-1 )/z * (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4])*n) ) -
		     2 * rho[2] * m


	          )
}

drhonn <- function(rho, parms, z = 4) { 
	  with(parms,
2 * (rho[5] - rho[3] - rho[1]) *
		     ( del*rho[5] + (1-del)/z + (1-del) * (z-1)/z * (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4]) )*
		     (b- cn* rho[5] - cpn* rho[4] ) - 2 * rho[3] * m  
	       )
}

drhop <- function(rho, parms, z = 4) { 
	  with(parms,                             
	(1- rho[5] - rho[4]) * (del * rho[4] +  (1-del) * (rho[4] - rho[2] - rho[1])/(1- rho[5] - rho[4])) *
		     (b- cp* rho[4] - cnp* rho[5] - g*( 1 - (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4])*n) ) - rho[4] * m 
		             )
}

drhon <- function(rho, parms, z = 4) { 
	  with(parms,
(1- rho[4] - rho[5]) * ( del*rho[5] + (1-del) * (rho[5] - rho[3] - rho[1])/(1- rho[5] - rho[4]) ) *
		     (b- cn* rho[5] - cpn* rho[4] ) - rho[5] * m  
		           )
}
# }}}

# Système d'équations différentielles# {{{
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

# Parameters definitions# {{{
rhoini <- 0.1

rho <- c(rhonp = rhoini*0.4, # rho[1]
	 rhopp = rhoini*0.4 , # rho[2]
	 rhonn = (1-rhoini)*0.6, # rho[3]
	 rhop  = rhoini,  # rho[4]
	 rhon  = (1-rhoini)*0.9   # rho[5]
	 )

# Parallel parameters setting
	 modelparms <- list(del = 0.1,
		   b   = seq(0.8, 0.4, -0.1),
		   cn  = seq(0.05, 0.2, 0.05),
		   cp  = seq(0.05, 0.2, 0.05),
		   cnp = seq(0.05, 0.2, 0.05),
		   cpn = seq(0.2, 0.4, 0.05),
		   g   = seq(0.0, 0.2, 0.05),
		   n   = seq(0, 1, 1),
		   m   = 0.2
	  )

iterations <- expand.grid(modelparms)  #Etabli toutes les combinaisons de paramètre possibles. 
iterations <- cbind(ID = 1:dim(iterations)[1],
                    iterations)
str(iterations)
# }}}


######################### output saving

dir.create("test_pairs")
setwd(paste(getwd(),"/test_pairs",sep=""))
######################### starting parallel backend

library(foreach)
library(doSNOW)
library(deSolve) 

workerlist <- rep("localhost", times = 23)
cl         <- makeSOCKcluster(workerlist)
registerDoSNOW(cl)

# Foreach loop# {{{
output <- foreach(iteration = iterations$ID,
		  .combine  = rbind,
		  .packages = "deSolve") %dopar% {

	parsTemp <- as.list(iterations[iterations$ID == iteration,]) #Take parameters of the iteration.
	print(parsTemp)

	# Run pair-approx simulation
	runmodel <- ode(y = rho, func = odesys, times = seq(1, 3000, 30), parms = parsTemp)

	out       <- as.data.frame(runmodel)
        out       <- round(out ,3)
        out$rho0  <- 1 - out$rhop - out$rhon
        out$rhop0 <- with(out, rhop-rhonp-rhopp)
        out$rhon0 <- with(out, rhon-rhonp-rhonn)
        out$rho00 <- with(out, rho0-rhon0-rhop0)
        out$Cnp   <- with(out, (rhonp/rhon)/rhop)

	result <- list(parsTemp, out)
# 	bien définir le output ci-dessous
	output <- data.frame(
			     ID   = iteration,
			     del  = parsTemp$del,
			     m    = parsTemp$m,
			     n    = parsTemp$n,
			     b    = parsTemp$b,
			     cn   = parsTemp$cn,
			     cp   = parsTemp$cp,
			     cnp  = parsTemp$cnp,
			     cpn  = parsTemp$cpn,
			     g    = parsTemp$g,
			     rhop = out$rhop[nrow(out)],
			     rhon = out$rhon[nrow(out)],
			     Cnp  = out$Cnp[nrow(out)]
			    ) 
			     
	save(result, file = paste("result-",iteration,".rdata", sep = ""))
	gc() 
	return(output) 
}
# }}}
write.csv(output, file = "output.csv",col.names =colnames(output))

stopCluster(cl)


