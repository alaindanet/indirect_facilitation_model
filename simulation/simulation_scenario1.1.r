##################################
#
# Author: Alain Danet
#
# From Florian Schneider
#
# Date: 25.04.2014
#
##################################



rm(list=ls())

########################################################################################

count <- function(x, neighbor) {
  x_logical_with_border <- (x$cells %in% neighbor)[x_with_border]
  rowSums(sapply(interact,function(x){x_logical_with_border[x_to_evaluate+x]}))
}


######################### defining parameter set

parameters = list(
  m1 = 0.2,   # intrinsic mortality
  #b1 = 0.62,  # beta*eps
  d = 0.1,		# degradation
  c_1 = 0.2, 		# beta*g
  del1 = 0.1, # seeds dispersed; (1-del) seeds on nearest neighbourhood
  r = 0.0, 	# regeneration rate
  f = 0.9 		# local facilitation
)

# Add parameters for second species

parameters <- c(
  parameters,
  starting = 0.5,
  m2 = as.numeric(parameters["m1"]),
  c_2 = as.numeric(parameters["c_1"]),
  del2 = as.numeric(parameters["del1"]),
  p=0.0,
  com="poly"
)

First_ID = 1

# Parallel parameters setting

parallel_param <- list(b1 = seq(0.7,0.5,-0.1),
                       g = seq(0.0,0.1,0.05),
                       n = seq(0,1,1),
                       c_21 = seq(0.05,0.2,0.05),
                       c_12 = seq(0.05,0.2,0.05)
)
iterations <- expand.grid(parallel_param)

iterations <- cbind(ID = 1:dim(iterations)[1],
                    iterations,
                    b2 = iterations$b1,
                    parameters)
str(iterations)

min_replicates = 10
max_tries = 50

# specify lattice
width = 50
height = 50

# initial cell states
states = c("+1","+2","0","-")# +1 will be the nurse, +2 the protege

color <- c("green","black","grey80", "white")


# time and resolution of simulation
#timesteps = 1000
delta = 1/5
t_min = 500
t_max = 2500
t_eval <- 200


################ map objects

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
neighbours_in_I <- which(is.finite(abs(I)/abs(I)), arr.in = TRUE)
# coordinates relative to the evaluated cell (=  which(is.na(I) ) 
relrow <- neighbours_in_I[,1]-which(is.na(I), arr.ind = TRUE)[1]
relcol <- neighbours_in_I[,2]-which(is.na(I), arr.ind = TRUE)[2]

# relative position of the four direct neighbours of a cell
interact <- (relrow * dim(X)[2] + relcol)

######################### output saving

dir.create("test_parallel")
setwd("/home/alain/test_parallel")
######################### starting parallel backend

library(foreach)
library(doSNOW)

workerlist <- rep("localhost", times = 10)

cl <- makeSOCKcluster(workerlist)

registerDoSNOW(cl)
# run parameters (in parallel)

foreach(iteration = iterations$ID, .combine = rbind) %dopar% { 
  
  
  set.seed(iterations$seed[iteration])
  
  # collecting final lattices for each run in j
  collect <- list()
  collect$lattice <- list()
  collect$out <- data.frame(
    replicate = NA,        
    
    rho_plus_ini = NA,
    rho_plus = NA, 
    rho_plus_sd = NA,
    rho_plus_fin = 0,
    
    rho_nurse_ini = NA,
    rho_nurse = NA, 
    rho_nurse_sd = NA,
    rho_nurse_fin = 0,
    
    rho_protege_ini = NA,
    rho_protege = NA, 
    rho_protege_sd = NA,
    rho_protege_fin = 0,
    
    q_1_1_ini = NA,
    q_1_1 = NA,
    q_1_1_sd = NA,
    q_1_1_fin = NA,
    
    q_2_2_ini = NA,
    q_2_2 = NA,
    q_2_2_sd = NA,
    q_2_2_fin = NA,
    
    q_1_2_ini = NA,
    q_1_2 = NA,
    q_1_2_sd = NA,
    q_1_2_fin = NA,
    
    q_2_1_ini = NA,
    q_2_1 = NA,
    q_2_1_sd = NA,
    q_2_1_fin = NA,
    
    clus_1_1 = NA,
    clus_1_1_sd = NA,
    
    clus_2_2 = NA,
    clus_2_2_sd = NA,
    
    clus_1_2 = NA,
    clus_1_2_sd = NA,
    
    clus_2_1 = NA,
    clus_2_1_sd = NA,
    
    stability_total = 0,
    stability_nurse = 0,
    stability_protege = 0,
    runtime = NA,
    starting = NA
  )[-1,]
  
  
  # run replicates
  
  j = 0
  n = 0
  
  while(n+1 <= min_replicates & j+1 <= max_tries) {
    
    
    parms_temp <- as.list(subset(iterations, ID == iteration))
    
    ###### initialise grid
    # how many initial plants
    init_plant <- as.integer(width*height*parms_temp$starting)
    
    # vector of empty cells in state "0"
    cells <- factor(rep("0", times = width*height), levels = states) # 1:length(states))
    # replace init_plant cells, randomly drawn, with "+"
    
    if(parms_temp$com == "mono_p"){
      cells[sample(1:(width*height), size = init_plant, replace = FALSE)] <- "+2"
      empty <- which(! cells %in% "+2")
      # select which will be degraded? fixed to 50% of non occupied cells. 
      init_degraded <- sample(empty, size = length(empty)/2) 
      # replace cell state. 
      cells[init_degraded] <- "-"
    }
      if(parms_temp$com == "mono_n"){
        cells[sample(1:(width*height), size = init_plant, replace = FALSE)] <- "+1"
        empty <- which(! cells %in% "+1")
        # select which will be degraded? fixed to 50% of non occupied cells. 
        init_degraded <- sample(empty, size = length(empty)/2) 
        # replace cell state. 
        cells[init_degraded] <- "-"
    }
    if(parms_temp$com == "poly"){
      cells[sample(1:(width*height), size =  init_plant/2, replace = FALSE)] <- "+1"
      cells[sample(which(cells != "+1"), size = init_plant/2, replace = FALSE)] <- "+2"
      
      # which cells are still empty?
      empty <- which(! cells %in% c("+1","+2"))
      # select which will be degraded? fixed to 50% of non occupied cells. 
      init_degraded <- sample(empty, size = length(empty)/2) 
      # replace cell state. 
      cells[init_degraded] <- "-"
    }
    
    initial <- list(  
      dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
      cells = cells#contains a random row-wise, factorial vector to fill the grid 
    )
    levels(initial$cells) <- states  #assign cell states 
    class(initial) <- c("list","landscape") # set class of object (required for plotting)
    parms_temp$rho_plus <- sum( sum(initial$cells == "+1"),sum(initial$cells == "+2") )/(width*height)    
    
    #### initialise result object
    result <- list()
    result$time <- seq(0, t_min, delta)
    
    
    result$rho_plus <- vector("numeric", length = length(result$time))
    result$rho_plus[1] <- parms_temp$rho_plus
    
    result$rho_nurse <- vector("numeric", length = length(result$time))
    result$rho_nurse[1] <- sum(initial$cells == "+1")/(width*height)
    
    result$rho_protege <- vector("numeric", length = length(result$time))
    result$rho_protege[1] <- sum(initial$cells == "+2")/(width*height)
    
    
    result$q_1_1  <- vector("numeric", length = length(result$time))  		
    result$q_1_1[1]  <- mean(subset(count(initial, "+1")/4, initial$cells == "+1") )
    
    result$q_2_2  <- vector("numeric", length = length(result$time))  		
    result$q_2_2[1]  <- mean(subset(count(initial, "+2")/4, initial$cells == "+2") )
    
    result$q_1_2  <- vector("numeric", length = length(result$time))    	
    result$q_1_2[1]  <- mean(subset(count(initial, "+1")/4, initial$cells == "+2") )
    
    result$q_2_1  <- vector("numeric", length = length(result$time))    	
    result$q_2_1[1]  <- mean(subset(count(initial, "+2")/4, initial$cells == "+1") )
    
    
    x_old <- initial
    
    stability_total <- 1
    stability_nurse <- 1
    stability_protege <- 1
    
    i = 1
    
    # starting iterations
    while(stability_total > 0.0001 | stability_nurse > 0.05 | stability_protege > 0.05 & i <= t_max/delta & result$rho_plus[i] > 0) {
      
      i <- i +1
      x_new <- x_old 		# copy x_old into an object x_new to allocate memory	
      
      # model specific part:
      # 1 - setting time-step parameters
      
      # count local density of occupied fields for each cell: 
      parms_temp$Q_plus1 <- count(x_old, "+1")/4
      parms_temp$Q_plus2 <- count(x_old, "+2")/4
      parms_temp$rho_nurse <- result$rho_nurse[i-1]#sum(x_old$cells == "+1")/(width*height)
      parms_temp$rho_protege <- result$rho_protege[i-1]
      
      # 2 - drawing random numbers
      rnum <- runif(width*height) # one random number between 0 and 1 for each cell
      
      # 4 - applying the rules to fill the cells in x_new
      
      # calculate recolonisation rates of all cells
      recolonisation1 <- with(parms_temp, (del1*rho_nurse+(1-del1)*Q_plus1)*(b1-c_1*rho_nurse-c_21*rho_protege-g*Q_plus2*p)*delta)
      recolonisation2 <- with(parms_temp, (del2*rho_protege+(1-del2)*Q_plus2)*(b2-c_2*rho_protege-c_12*rho_nurse-g*(1-Q_plus1*n))*delta)
      # calculate death rates                              
      death1 <- with(parms_temp, m1*delta)
      death2 <- with(parms_temp, m2*delta)
      
      # calculate regeneration rate and degradation rate
      regeneration <- with(parms_temp, (r + f*(Q_plus1+Q_plus2)*delta))
      degradation <- with(parms_temp, (d *delta))
      
      # check for sum of probabilities to be inferior 1 and superior 0
      if(any(c(recolonisation1+degradation,
               recolonisation2+degradation,  recolonisation1 + recolonisation2 + degradation,
               death1, death2, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in run",
                                                                  iteration, "time step", i, "! decrease delta!!!")) 
      if(any(c(recolonisation1, recolonisation2,
               degradation, regeneration,
               death1, death2) < 0)) warning(paste("a set probability falls below 0 in run",
                                                   iteration, "in time step", i, "! balance parameters!!!")) 
      
      # apply rules 
      
      x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation1)] <- "+1"
      x_new$cells[which(x_old$cells == "0" & rnum > recolonisation1 & rnum <= recolonisation1+recolonisation2)] <- "+2"
      x_new$cells[which(x_old$cells == "0" & rnum > recolonisation1 + recolonisation2 & rnum <= recolonisation1 + recolonisation2 + degradation)] <- "-"
      x_new$cells[which(x_old$cells == "+1" & rnum <= death1)] <- "0"
      x_new$cells[which(x_old$cells == "+2" & rnum <= death2)] <- "0"
      x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"
      
      # 5 saving state of the new grid		
      
       
      result$rho_plus[i] <- sum( sum(x_new$cells == "+1"), sum(x_new$cells == "+2") )/(width*height)
      
      result$rho_nurse[i] <- sum(x_new$cells == "+1")/(width*height)
      
      result$rho_protege[i] <- sum(x_new$cells == "+2")/(width*height)
      
      result$q_1_1[i]  <- mean(subset(count(x_new, "+1")/4, x_new$cells == "+1") )

      result$q_2_2[i]  <- mean(subset(count(x_new, "+2")/4, x_new$cells == "+2") )
      
      result$q_1_2[i]  <- mean(subset(count(x_new, "+1")/4, x_new$cells == "+2") )

      result$q_2_1[i]  <- mean(subset(count(x_new, "+2")/4, x_new$cells == "+1") )
      
      
      x_old <- x_new
      
      
      if(i > t_min/delta+1 ) { #| result$rho_plus[i] == 0
     
        
        if( result$rho_plus[i] > 0) {
          t_1 <- (i-2*t_eval/delta):(i-t_eval/delta)-1
          t_2 <- (i-t_eval/delta):(i)
          stability_total <- (abs(mean(result$rho_plus[t_1]) - mean(result$rho_plus[t_2])))/(mean(result$rho_plus[t_1]))
          
          if( result$rho_nurse[i] > 0){
          stability_nurse <- (abs(mean(result$rho_nurse[t_1]) - mean(result$rho_nurse[t_2])))/(mean(result$rho_nurse[t_1]))
          } else {stability_nurse = 0}
          
          if( result$rho_protege[i] > 0){
          stability_protege <- (abs(mean(result$rho_protege[t_1]) - mean(result$rho_protege[t_2])))/(mean(result$rho_protege[t_1]))
          } else {stability_protege = 0}
        
          } else {
          stability_total = 0
          stability_nurse = 0
          stability_protege = 0
        }
        
        result$time[i] <- i*delta          
        
      }

    } # end of simulation run (over i)
    
    j = (j + 1)
    
    if(stability_total <= 0.0001 & stability_nurse <= 0.05 & stability_protege <= 0.05) { 
      n = (n + 1) 
      
      collect$lattice[[n]] <- x_new
      
      t_fin = result$time[length(result$time)]
      i_stable <- (length(result$time)-(2*t_eval/delta)):length(result$time)
      
      # switch for output at extinction
      if(result$rho_plus[t_fin] > 0) {
        result$out <- data.frame(
          
          replicate = j,
          
          rho_plus_ini = result$rho_plus[1],
          rho_plus = mean(result$rho_plus[i_stable], na.rm=T), 
          rho_plus_sd = sd(result$rho_plus[i_stable], na.rm=T),
          rho_plus_fin = result$rho_plus[length(result$time)],
          
          rho_nurse_ini = result$rho_nurse[1],
          rho_nurse = mean(result$rho_nurse[i_stable], na.rm=T), 
          rho_nurse_sd = sd(result$rho_nurse[i_stable], na.rm=T),
          rho_nurse_fin = result$rho_nurse[length(result$time)],
          
          rho_protege_ini = result$rho_protege[1],
          rho_protege = mean(result$rho_protege[i_stable], na.rm=T), 
          rho_protege_sd = sd(result$rho_protege[i_stable], na.rm=T),
          rho_protege_fin = result$rho_protege[length(result$time)],
          
          q_1_1_ini = result$q_1_1[1],
          q_1_1 = mean(result$q_1_1[i_stable], na.rm=T),
          q_1_1_sd = sd(result$q_1_1[i_stable], na.rm=T),
          q_1_1_fin = result$q_1_1[length(result$time)],
          
          q_2_2_ini = result$q_2_2[1],
          q_2_2 = mean(result$q_2_2[i_stable], na.rm=T),
          q_2_2_sd = sd(result$q_2_2[i_stable], na.rm=T),
          q_2_2_fin = result$q_2_2[length(result$time)],
          
          q_1_2_ini = result$q_1_2[1],
          q_1_2 = mean(result$q_1_2[i_stable], na.rm=T),
          q_1_2_sd = sd(result$q_1_2[i_stable], na.rm=T),
          q_1_2_fin = result$q_1_2[length(result$time)],
          
          q_2_1_ini = result$q_2_1[1],
          q_2_1 = mean(result$q_2_1[i_stable], na.rm=T),
          q_2_1_sd = sd(result$q_2_1[i_stable], na.rm=T),
          q_2_1_fin = result$q_2_1[length(result$time)],
          
          clus_1_1 = mean(result$q_1_1[i_stable]/result$rho_nurse[i_stable], na.rm=T),
          clus_1_1_sd = sd(result$q_1_1[i_stable]/result$rho_nurse[i_stable], na.rm=T),
          
          clus_2_2 = mean(result$q_2_2[i_stable]/result$rho_protege[i_stable], na.rm=T),
          clus_2_2_sd = sd(result$q_2_2[i_stable]/result$rho_protege[i_stable], na.rm=T),
          
          clus_1_2 = mean(result$q_1_2[i_stable]/result$rho_nurse[i_stable], na.rm=T),
          clus_1_2_sd = sd(result$q_1_2[i_stable]/result$rho_nurse[i_stable], na.rm=T),
          
          clus_2_1 = mean(result$q_2_1[i_stable]/result$rho_protege[i_stable], na.rm=T),
          clus_2_1_sd = sd(result$q_2_1[i_stable]/result$rho_protege[i_stable], na.rm=T),
          
          stability_total = stability_total,
          stability_nurse = stability_nurse,
          stability_protege = stability_protege,
          runtime = t_fin,
          starting = init_plant
        )
        
      } else {
        
        result$out <- data.frame(
          
          replicate = j,        
          
          rho_plus_ini = result$rho_plus[1],
          rho_plus = 0, 
          rho_plus_sd = 0,
          rho_plus_fin = 0,
          
          rho_nurse_ini = result$rho_nurse[1],
          rho_nurse = 0, 
          rho_nurse_sd = 0,
          rho_nurse_fin = 0,
          
          rho_protege_ini = result$rho_protege[1],
          rho_protege = 0, 
          rho_protege_sd = 0,
          rho_protege_fin = 0,
          
          q_1_1_ini = result$q_1_1[1],
          q_1_1 = NA,
          q_1_1_sd = NA,
          q_1_1_fin = NA,
          
          q_2_2_ini = result$q_2_2[1],
          q_2_2 = NA,
          q_2_2_sd = NA,
          q_2_2_fin = NA,
          
          q_1_2_ini = result$q_1_2[1],
          q_1_2 = NA,
          q_1_2_sd = NA,
          q_1_2_fin = NA,
          
          q_2_1_ini = result$q_2_1[1],
          q_2_1 = NA,
          q_2_1_sd = NA,
          q_2_1_fin = NA,
          
          clus_1_1 = NA,
          clus_1_1_sd = NA,
          
          clus_2_2 = NA,
          clus_2_2_sd = NA,
          
          clus_1_2 = NA,
          clus_1_2_sd = NA,
          
          clus_2_1 = NA,
          clus_2_1_sd = NA,
          
          
          stability_total = 0,
          stability_nurse = 0,
          stability_protege = 0,
          runtime = t_fin,
          starting = init_plant
        ) 
      }
      
      collect$out[n,] <- result$out
    }
   
    
  }#-> collect$out
  
  
  result <- list()
  
  # pool replicates to one row (means, sd)
  stable <- collect$out$stability_total < 0.0001 & collect$out$stability_nurse < 0.05 & collect$out$stability_protege < 0.05
  
  result$runs <- collect$out 
  result$grids <- collect$lattice
  
  
  result$out <- data.frame(
    ID = iteration,
    g = parms_temp$g, 
    b = parms_temp$b1, 
    m = parms_temp$m1,
    com = parms_temp$com,
    starting = parms_temp$starting,
    c_1 = parms_temp$c_1,
    c_2 = parms_temp$c_2,
    c_12 = parms_temp$c_12,
    c_21 = parms_temp$c_21,
    n = parms_temp$n,
    n_rep = length(which(stable)),
    runtime = mean(result$runs$runtime[stable], na.rm=T), 
    runtime_sd = sd(result$runs$runtime[stable], na.rm=T), 
    
    rho_plus_ini = mean(result$runs$rho_plus_ini[stable], na.rm=T),
    rho_plus = mean(result$runs$rho_plus[stable], na.rm=T),
    rho_plus_sd = mean(result$runs$rho_plus_sd[stable], na.rm=T),
    rho_plus_sd_inter_replicates = sd(result$runs$rho_plus[stable], na.rm=T),
    
    rho_nurse_ini = mean(result$runs$rho_nurse_ini[stable], na.rm=T),
    rho_nurse = mean(result$runs$rho_nurse[stable], na.rm=T), 
    rho_nurse_sd = mean(result$runs$rho_nurse_ini[stable], na.rm=T),
    rho_nurse_sd_inter_replicates = sd(result$runs$rho_nurse[stable], na.rm=T),
    
    rho_protege_ini = mean(result$runs$rho_protege_ini[stable], na.rm=T),
    rho_protege = mean(result$runs$rho_protege[stable], na.rm=T), 
    rho_protege_sd = mean(result$runs$rho_protege_ini[stable], na.rm=T),
    rho_nurse_sd_inter_replicates = sd(result$runs$rho_protege[stable], na.rm=T),
    
    q_1_1 =  mean(result$runs$q_1_1[stable], na.rm=T),
    q_1_1_sd = mean(result$runs$q_1_1_sd[stable], na.rm=T),
    q_1_1_inter_replicates = sd(result$runs$q_1_1[stable], na.rm=T),
    
    q_2_2 =  mean(result$runs$q_2_2[stable], na.rm=T),
    q_2_2_sd = mean(result$runs$q_2_2_sd[stable], na.rm=T),
    q_2_2_inter_replicates = sd(result$runs$q_2_2[stable], na.rm=T),
    
    q_1_2 =  mean(result$runs$q_1_2[stable], na.rm=T),
    q_1_2_sd = mean(result$runs$q_1_2_sd[stable], na.rm=T),
    q_1_2_inter_replicates = sd(result$runs$q_1_2[stable], na.rm=T),
    
    q_2_1 =  mean(result$runs$q_2_1[stable], na.rm=T),
    q_2_1_sd = mean(result$runs$q_2_1_sd[stable], na.rm=T),
    q_2_1_inter_replicates = sd(result$runs$q_2_1[stable], na.rm=T),
    
    clus_1_1 = mean(result$runs$clus_1_1[stable], na.rm=T),
    clus_1_1_sd = sd(result$runs$clus_1_1_sd[stable], na.rm=T),
    clus_1_1_inter_replicates = sd(result$runs$clus_1_1[stable], na.rm=T),
    
    clus_2_2 = mean(result$runs$clus_2_2[stable], na.rm=T),
    clus_2_2_sd = sd(result$runs$clus_2_2_sd[stable], na.rm=T),
    clus_2_2_inter_replicates = sd(result$runs$clus_2_2[stable], na.rm=T),
    
    clus_1_2 = mean(result$runs$clus_1_2[stable], na.rm=T),
    clus_1_2_sd = sd(result$runs$clus_1_2_sd[stable], na.rm=T),
    clus_1_2_inter_replicates = sd(result$runs$clus_1_2[stable], na.rm=T),
    
    clus_2_1 = mean(result$runs$clus_2_1[stable], na.rm=T),
    clus_2_1_sd = sd(result$runs$clus_2_1_sd[stable], na.rm=T),
    clus_2_1_inter_replicates = sd(result$runs$clus_2_1[stable], na.rm=T),
    
    runtime = t_fin
  )
  
  
  
  
  save(result, file = paste("result_scenario1.1_", iteration,".rdata", sep = ""))
  
  gc() 
  return(result$out)
  
}-> output

write.csv(output, file = paste("output_scenario1.1.csv", sep = "") )

stopCluster(cl)
