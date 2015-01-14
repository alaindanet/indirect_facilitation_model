#################################################
#
#		model: facilitation model (Kefi et al 2007)
#
#
#author: Flo
#date: 02.04.2013
#
#
#	Tasks for Flo: 
#		identify bottlenecks for speedy and parallel computing
#		implement test for stable equilibrium
#		
#   Tasks for Marina:
# 		adapt for two-species model (implement rules, find parameters)
# 		adapt output
# 
#  
#################################################


rm(list=ls())

#
#setwd("sftp://alain@162.38.184.118/home/alain/test")

# plotting function for objects of class "landscape".
plot.landscape <- function(x, grid = FALSE, axis = FALSE, cols = "auto", add = FALSE, ani = FALSE, ...) {
  lvls <- levels(x$cells) 
  nlev <- length(lvls)
  if(cols[1] == "auto") cols = greyscale(nlev) # default color value
  
  if(ani & Sys.info()[['sysname']] == "Windows") adj = -0.5 else adj = 0 #this adjustment constant is added when producing a pixel accurate png or gif file, requires TRUE when the function is used to plot animated figures. 
  
  if(!add) plot(NA,NA, xlim = c(0.5+adj, x$dim[1]+0.5+adj), ylim = c( x$dim[2]+0.5+adj, 0+0.5+adj), bty = c("n", "o")[grid+1], xaxs = "i", yaxs = "i",xlab = "", ylab = "", xaxt = "n", yaxt = "n", ... ) 
  
  if(axis && !add) axis(3) 
  if(axis && !add) axis(2)
  
  if(grid) border = "grey80" else border = cols[as.numeric(x$cells)]
  
  rect(rep(1:x$dim[1], times = x$dim[2])-.5, rep(1:x$dim[2], each = x$dim[1])-.5, rep(1:x$dim[1], times = x$dim[2])+.5, rep(1:x$dim[2], each = x$dim[1])+.5, col = cols[as.numeric(x$cells)], border = border)
  
  if(grid) box()
 }


count  <- function(x, neighbor) {
  x_logical_with_border <- (x$cells %in% neighbor)[x_with_border]
  rowSums(sapply(interact,function(x){x_logical_with_border[x_to_evaluate+x]}))
}




# specify lattice
width = 50
height = 50



# time and resolution of simulation
timesteps = 500
delta = 1/5


# derive helper vectors for counting: 
# transformation vector for evaluation at the border of the grid
	# set evaluation matrix 
	X <- matrix(as.integer(1:(width*height)), ncol = width, byrow =TRUE)
	# setting the border of the evaluation matrix X
	X <- cbind(X[,width], X, X[,1] )  
	X <- rbind(X[height,], X, X[1,] ) 
	# transformation vector which adds the border to the lattice:
	x_with_border <- as.integer(t(X))
	
	# from the matrix X (with copied border cells), which cells are the actual cells (reverse transformation vector of the previous lines) 
	x_to_evaluate <- sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]	)		

# defining the neighborhood which is to be evaluated	
	# set interaction matrix
	I <- matrix(c(0,1,0,1,NA,1,0,1,0), ncol = 3, byrow = TRUE)	
	# coordinates of neighbours in Interaction matrix I: 
	neighbours_in_I <- which(is.finite(abs(I)/abs(I)), arr.ind = TRUE)
	# coordinates relative to the evaluated cell (=  which(is.na(I) ) 
	relrow <- neighbours_in_I[,1]-which(is.na(I), arr.ind = TRUE)[1]
	relcol <- neighbours_in_I[,2]-which(is.na(I), arr.ind = TRUE)[2]
	
	# relative position of the four direct neighbours of a cell
	interact <- (relrow * dim(X)[2] + relcol)

# initial cell states
states = c("+1","+2","0","-")# +1 will be the nurse, +2 the protege
prob = c(0.25,0.25,0.25,0.25)	

# sampling the initial random grid into a list object
initial <- list(  
	dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
	cells = sample(factor(1:length(states)), width*height, replace = T, prob = prob ) #second element contains a random row-wise, factorial vector to fill the grid 
	)
levels(initial$cells) <- states  #assign cell states 
class(initial) <- c("list","landscape") # set class of object (required for plotting)

color <- c("green1","black","grey80", "white") # define colors for the cell state levels; green: nurse; black:protege

#plot the initial state
par(mar = c(0,0,0,0), mfrow = c(1,1))
#x11()
plot(initial, cols = color)

# defining parameter set
parameters = list(
	m1 = 0.2, m2 = 0.2,  # intrinsic mortality
	b1 = 0.8,  # beta*eps
	d = 0.1,		# degradation
	c_1 = 0.1, c_2 = 0.1,		# beta*g
	c_12 = 0.05, c_21 = 0.2, 
	del1 = 0.1, # seeds dispersed; (1-del) seeds on nearest neighbourhood
	r = 0.0, 	# regeneration rate
	f = 0.9,  # local fascilitation
	g = 0.1,
	n = 0, # protege protection against herbivory (associational résistance)
	p = 0,
  g2 = 0,# grazing
  al1 = 0,
  al2 = 0,
  n1 = 0,
  n2 = 0
  )
# Add parameters for second species
parameters <- c(
  parameters,
  b2 = as.numeric(parameters["b1"]),
  del2 = as.numeric(parameters["del1"])
  )

  

# initialising result list 
result <- list()   # create an output list and
			result$time <- seq(0, timesteps, delta) # write simulated timesteps
			
			#save global densities of state i: rho_i
			result$rho  <- list()			
			#preallocate memory for the saving
			for(i in 1:length(states)) {
				result$rho[[i]] <- numeric(length = timesteps/delta)
				result$rho[[i]][1]  <- sum(initial$cells == states[i])/(width*height) # write initial rho 
			}
			
			#save local densities of state i given that a cell is in state j: q_i_j
			result$q_1  <- list()		
			#preallocate memory for the saving	
			for(j in c(1:2)) {
			result$q_1[[j]] <- numeric(length = timesteps/delta)
			result$q_1[[j]][1]  <- mean(rowSums( sapply(interact,
                                                  function(k) (initial$cells == states[j])[x_with_border][x_to_evaluate[initial$cells == states[1]]+k]) )) /4# write initial rho 
			}
			
      result$q_2 <- list()
      
      for(j in 1:2) {
      result$q_2[[j]] <- numeric(length = timesteps/delta)
      result$q_2[[j]][1]  <- mean(rowSums( sapply(interact,
                                                  function(k) (initial$cells == states[j])[x_with_border][x_to_evaluate[initial$cells == states[2]]+k]) )) /4# write initial rho 
      }

			# save only relevant q
			
			result$timeseries <- list() # create a subordinate list and
			for(i in seq(1,timesteps+1, 1)) result$timeseries[[i]] <- initial	# allocate memory for each timeseries object
			

#  initialising first iteration 
	x_old <- initial    #old landscape for first timestep   
	parms_temp <- parameters # copying parms for the simulation and multiplying
	
for(i in seq_along(result$time)) {    #calculation loop

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho1 <- result$rho[[1]][i-1]   # get old rho plus for this timestep 
		parms_temp$rho2 <- result$rho[[2]][i-1]
	 # count local density of occupied fields for each cell: 
 		parms_temp$Q_plus1 <- count(x_old, "+1")/4
    parms_temp$Q_plus2 <- count(x_old, "+2")/4
		# This is a vectorised evaluation: first transform the original grid x_old, into a grid with border. Here check for all cells of the original grid (x_to_evaluate) the content of all four neighbouring cells, whose relative position is defined in interact. The number of cells with content "+" are count and divided by 4.
	# count local density of degraded fields for each cell
		#parms_temp$Q_minus <- rowSums(sapply(interact, function(j)  (x_old$cells == "-")[x_with_border][x_to_evaluate+j]))/4 
		
# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
		recolonisation1 <- with(parms_temp, (del1*rho1+(1-del1)*Q_plus1)*(b1-c_1*(rho1 + rho2) - c_21*Q_plus1 - g*Q_plus2*p)*delta)
    recolonisation2 <- with(parms_temp, (del2*rho2+(1-del2)*Q_plus2)*(b2-c_2*(rho1 + rho2) - c_12*Q_plus2 - g*(1-Q_plus1*n))*delta)
		degradation <- with(parms_temp, (d *delta))
		death1 <- with(parms_temp, (m1 + g2*(al1 - n1*Q_plus1))*delta)
    death2 <- with(parms_temp, (m2 + g2*(al2 - n2*Q_plus1))*delta)
		regeneration <- with(parms_temp, (r + f*(Q_plus1+Q_plus2))*delta)
		
		# maybe find a way to calculate only relevant probs
		
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation1)] <- "+1"
    x_new$cells[which(x_old$cells == "0" & rnum > recolonisation1 & rnum <= recolonisation1+recolonisation2)] <- "+2"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation1 + recolonisation2 & rnum <= recolonisation1 + recolonisation2 + degradation)] <- "-"
		x_new$cells[which(x_old$cells == "+1" & rnum <= death1)] <- "0"
    x_new$cells[which(x_old$cells == "+2" & rnum <= death2)] <- "0"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"

	# put warning if 'recolonisation1 + recolonisation2 + degradation' > 1
if(any(c(recolonisation1 + recolonisation2 + degradation) > 1 )) warning(paste("a set probability is exceeding 1 in run",
                                                                               "time step", i, "! decrease delta!!!")) 
if(any(c(recolonisation1, recolonisation2, degradation, regeneration) < 0)) warning(paste("a set probability falls below 0 in run",
                                                                               "in time step", i, "! balance parameters!!!"))

# 5 saving state of the new grid		
		for(j in 1:2){
			result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
			result$q_1[[j]][i]  <-  mean(subset(count(x_new, states[j]), x_new$cells == "+1") /4)
			result$q_2[[j]][i]  <-  mean(subset(count(x_new, states[j]), x_new$cells == "+2") /4)
			}

		
		#for(j in 1:length(states)) {
		#result$q_[[j]][i]  <-  mean(subset(count(x_new, states[j], "+"), x_new$cells == states[j]) /4)
		#}
		
		# activate to save each single timeseries step
		#if(result$time[i] %in% seq(1,timesteps+1, 1)) result$timeseries[[ result$time[i] ]] <- x_new  #the whole grid is saved to timeseries
		
		# activate for plotting during simulation (very slow !!!)
		#if(result$time[i] %in% seq(1,timesteps+1, 1)) plot(x_new, col = color, add = TRUE)
		x_old <- x_new 
#plot(x_new, cols = color, add= TRUE)
		#gc()  #garbage collection
	} # end of simulation.


	
	
	
# the rest is graphical output

# FIGURE 1 	-- final states of the grid
par(mar = c(0,0,0,0))
	plot(x_new, cols = color)


# FIGURE 2 	-- global densities and local densities around cells in state "+" 
layout(matrix(1:2, ncol = 2), width = c(1,1))
	par(mar = c(4,4,2,1)+.1, xaxt = "s", yaxt ="s")

	plot(NA,NA, ylim = c(0, 1), xlim = range(result$time), type ="l", xlab = "Temps", ylab = expression(rho), xaxs ="i", yaxs = "i" )

		lines(result$time, result$rho[[1]])

		polygon(c(0,result$time,max(result$time)), c(0, result$rho[[1]]+result$rho[[2]]+result$rho[[3]]+result$rho[[4]], 0), col = "white")
		polygon(c(0,result$time,max(result$time)), c(0, result$rho[[1]]+result$rho[[2]]+result$rho[[3]], 0), col = "grey50")
		polygon(c(0,result$time,max(result$time)), c(0, result$rho[[1]]+result$rho[[2]], 0), col = "black")
polygon(c(0,result$time,max(result$time)), c(0, result$rho[[1]], 0), col = "green")

plot(result$rho[[2]] ~ result$time, ylim = c(0, 1), xlim = range(result$time),
     type ="l", lwd = 2, xlab = "Temps", ylab = expression(rho), xaxs ="i", yaxs = "i")
lines(result$rho[[1]]~ result$time,  lwd = 2, col="green1")

plot( I(result$q_1[[2]]/result$rho[[2]]) ~ result$time, ylim = c(0, 1.5), xlim = range(result$time),
      type ="l", lwd = 2, xlab = "Temps", ylab = expression(q["+1+2"]), xaxs ="i", yaxs = "i")
lines(I(result$rho[[2]]) ~ result$time, lwd = 2, col="red")
	plot(x_new, col = color)



# FIGURE 3 -- animated gif
library(animation)
if(Sys.info()[['sysname']] == "Linux") X11.options(antialias = "none") #for Linux Systems to enable pixel-wise plotting in (animated) gif-files. 
if(Sys.info()[['sysname']] == "Windows") windows.options(antialias = "none") #for Windows Systems to enable pixel-wise plotting in (animated) gif-files. 

saveGIF( 
	for(i in seq(1, length(result$timeseries), 1)) {
	    par(mar = c(0,0,0,0))
		plot(result$timeseries[[i]], grid = FALSE, cols = color, ani = TRUE)
	}
, movie.name = "facilitation_test.gif", img.name = "grid", convert = "convert", interval = 0.01/1,
    cmd.fun = system, clean = TRUE, ani.width = width, ani.height = height, outdir = getwd())
