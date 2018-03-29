#'#'###################################
#'  Model CA with local competition  #'
#'####################################'


rm(list=ls())

########################################################################################

count <- function(x, neighbor) {
  x.logical.with.border <- (x$cells %in% neighbor)[x.with.border]
  rowSums(sapply(interact,function(x){x.logical.with.border[x.to.evaluate+x]}))
}


# defining parameter set# {{{

parameters = list(m = 0.1,   # intrinsic mortality
		  d = 0.1,	# degradation
		  c1 = seq(-0.2,0.2,0.05), 		# beta*g
		  c2 = seq(-0.2,0.2,0.05),
		  del = 0.1, # seeds dispersed; (1-del) seeds on nearest neighbourhood
		  r = 0.01, 	# regeneration rate
		  f = 0.9, 		# local facilitation
		  starting = 0.5,
		  p=0.0,
		  com="poly",
		  b = seq(0.8,0.3,-0.1),
                  g = seq(0.0,0.3,0.05),
                  n = seq(0,1,1),
                  c21 = seq(-0.2,0.2,0.05),
                  c12 = seq(-0.2,0.2,0.05)
		  )

First.ID = 1



# Parallel parameters setting
iterations <- expand.grid(parameters)
iterations <- cbind(ID = 1:dim(iterations)[1],
                    iterations)
iterations  <- iterations[-which(iterations$c12 > iterations$c21),]

str(iterations)

min.replicates = 10
max.tries = 50

# specify lattice
width = 50
height = 50

# initial cell states
states = c("+1","+2","0","-")# +1 will be the nurse, +2 the protege

color <- c("green","black","grey80", "white")


# time and resolution of simulation
#timesteps = 1000
delta = 1/5
t.min = 500
t.max = 2500
t.eval <- 200
# }}}

# map objects# {{{

# derive helper vectors for counting: 
# transformation vector for evaluation at the border of the grid
# set evaluation matrix 
X <- matrix(as.integer(1:(width*height)), ncol = width, byrow =TRUE)
# setting the border of the evaluation matrix X
X <- cbind(X[,width], X, X[,1] )  
X <- rbind(X[height,], X, X[1,] ) 
# transformation vector which adds the border to the lattice:
x.with.border <- as.integer(t(X))

# from the matrix X (with copied border cells), which cells are the actual cells (reverse transformation vector of the previous lines) 
x.to.evaluate <- sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]	)		

# defining the neighborhood which is to be evaluated	
# set interaction matrix
I <- matrix(c(0,1,0,1,NA,1,0,1,0), ncol = 3, byrow = TRUE)	
# coordinates of neighbours in Interaction matrix I: 
neighbours.in.I <- which(is.finite(abs(I)/abs(I)), arr.in = TRUE)
# coordinates relative to the evaluated cell (=  which(is.na(I) ) 
relrow <- neighbours.in.I[,1]-which(is.na(I), arr.ind = TRUE)[1]
relcol <- neighbours.in.I[,2]-which(is.na(I), arr.ind = TRUE)[2]

# relative position of the four direct neighbours of a cell
interact <- (relrow * dim(X)[2] + relcol)
# }}}
######################### output saving

dir.create("CA17.05")
setwd(paste(getwd(),"/CA17.05/", sep="")) 
# starting parallel backend# {{{

library(foreach)
library(doSNOW)
#lecture du fichier contenant la liste des machines
args <- commandArgs(TRUE)
peFile=args[1]
liste <- Nodes=read.table(peFile, sep= "" ,header=F, stringsAsFactors=F)

#construction de la liste des slots pour le « cluster »

node <- Names=liste <- Nodes[,1]
nb <- slots=liste <- Nodes[,2]
workers=rep(node <- Names, nb <- slots)

nbworkers <- length(workers)

#We will Run in parallel mode (socket) with
cl <- makeSOCKcluster(workers)
registerDoSNOW(cl)

# run parameters (in parallel)

foreach(iteration = iterations$ID, .combine = rbind) %dopar% { 
  
  
  set.seed(iterations$seed[iteration])
  # }}}
  # collecting final lattices for each run in j# {{{
  collect <- list()
  collect$lattice <- list()
  collect$out <- data.frame(
    replicate = NA,        
    
    rho.plus.ini = NA,
    rho.plus = NA, 
    rho.plus.sd = NA,
    rho.plus.fin = 0,
    
    rho.nurse.ini = NA,
    rho.nurse = NA, 
    rho.nurse.sd = NA,
    rho.nurse.fin = 0,
    
    rho.protege.ini = NA,
    rho.protege = NA, 
    rho.protege.sd = NA,
    rho.protege.fin = 0,
    
    q.1.1.ini = NA,
    q.1.1 = NA,
    q.1.1.sd = NA,
    q.1.1.fin = NA,
    
    q.2.2.ini = NA,
    q.2.2 = NA,
    q.2.2.sd = NA,
    q.2.2.fin = NA,
    
    q.1.2.ini = NA,
    q.1.2 = NA,
    q.1.2.sd = NA,
    q.1.2.fin = NA,
    
    q.2.1.ini = NA,
    q.2.1 = NA,
    q.2.1.sd = NA,
    q.2.1.fin = NA,
    
    clus.1.1 = NA,
    clus.1.1.sd = NA,
    
    clus.2.2 = NA,
    clus.2.2.sd = NA,
    
    clus.1.2 = NA,
    clus.1.2.sd = NA,
    
    clus.2.1 = NA,
    clus.2.1.sd = NA,
    
    stability.total = 0,
    stability.nurse = 0,
    stability.protege = 0,
    runtime = NA,
    starting = NA
  )[-1,]# }}}
  
  
  # run replicates# {{{
  
  j = 0
  nrep = 0
  
  while(nrep+1 <= min.replicates & j+1 <= max.tries) {
    
    
    parms.temp <- as.list(subset(iterations, ID == iteration))
    
    ###### initialise grid
    # how many initial plants
    init.plant <- as.integer(width*height*parms.temp$starting)
    
    # vector of empty cells in state "0"
    cells <- factor(rep("0", times = width*height), levels = states) # 1:length(states))
    # replace init.plant cells, randomly drawn, with "+"
    
    if(parms.temp$com == "mono.p"){
      cells[sample(1:(width*height), size = init.plant, replace = FALSE)] <- "+2"
      empty <- which(! cells %in% "+2")
      # select which will be degraded? fixed to 50% of non occupied cells. 
      init.degraded <- sample(empty, size = length(empty)/2) 
      # replace cell state. 
      cells[init.degraded] <- "-"
    }
      if(parms.temp$com == "mono.n"){
        cells[sample(1:(width*height), size = init.plant, replace = FALSE)] <- "+1"
        empty <- which(! cells %in% "+1")
        # select which will be degraded? fixed to 50% of non occupied cells. 
        init.degraded <- sample(empty, size = length(empty)/2) 
        # replace cell state. 
        cells[init.degraded] <- "-"
    }
    if(parms.temp$com == "poly"){
      cells[sample(1:(width*height), size =  init.plant/2, replace = FALSE)] <- "+1"
      cells[sample(which(cells != "+1"), size = init.plant/2, replace = FALSE)] <- "+2"
      
      # which cells are still empty?
      empty <- which(! cells %in% c("+1","+2"))
      # select which will be degraded? fixed to 50% of non occupied cells. 
      init.degraded <- sample(empty, size = length(empty)/2) 
      # replace cell state. 
      cells[init.degraded] <- "-"
    }
    
    initial <- list(  
      dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
      cells = cells#contains a random row-wise, factorial vector to fill the grid 
    )
    levels(initial$cells) <- states  #assign cell states 
    class(initial) <- c("list","landscape") # set class of object (required for plotting)
    parms.temp$rho.plus <- sum( sum(initial$cells == "+1"),sum(initial$cells == "+2") )/(width*height)    # }}}
    
    #### initialise result object# {{{
    result <- list()
    result$time <- seq(0, t.min, delta)
    
    
    result$rho.plus <- vector("numeric", length = length(result$time))
    result$rho.plus[1] <- parms.temp$rho.plus
    
    result$rho.nurse <- vector("numeric", length = length(result$time))
    result$rho.nurse[1] <- sum(initial$cells == "+1")/(width*height)
    
    result$rho.protege <- vector("numeric", length = length(result$time))
    result$rho.protege[1] <- sum(initial$cells == "+2")/(width*height)
    
    
    result$q.1.1  <- vector("numeric", length = length(result$time))  		
    result$q.1.1[1]  <- mean(subset(count(initial, "+1")/4, initial$cells == "+1") )
    
    result$q.2.2  <- vector("numeric", length = length(result$time))  		
    result$q.2.2[1]  <- mean(subset(count(initial, "+2")/4, initial$cells == "+2") )
    
    result$q.1.2  <- vector("numeric", length = length(result$time))    	
    result$q.1.2[1]  <- mean(subset(count(initial, "+1")/4, initial$cells == "+2") )
    
    result$q.2.1  <- vector("numeric", length = length(result$time))    	
    result$q.2.1[1]  <- mean(subset(count(initial, "+2")/4, initial$cells == "+1") )
    
    
    x.old <- initial
    
    stability.total <- 1
    stability.nurse <- 1
    stability.protege <- 1
    
    i = 1
    # }}}
    # starting iterations# {{{
    while(stability.total > 0.0001 | stability.nurse > 0.05 | stability.protege > 0.05 & i <= t.max/delta & result$rho.plus[i] > 0) {
      
      i <- i +1
      x.new <- x.old 		# copy x.old into an object x.new to allocate memory	
      
      # model specific part:
      # 1 - setting time-step parameters
      
      # count local density of occupied fields for each cell: 
      parms.temp$Q.plus1 <- count(x.old, "+1")/4
      parms.temp$Q.plus2 <- count(x.old, "+2")/4
      parms.temp$rho.nurse <- result$rho.nurse[i-1]#sum(x.old$cells == "+1")/(width*height)
      parms.temp$rho.protege <- result$rho.protege[i-1]
      
      # 2 - drawing random numbers
      rnum <- runif(width*height) # one random number between 0 and 1 for each cell
      
      # 4 - applying the rules to fill the cells in x.new
      
      #if(flag) { # if density is unequal 0, then
      #recolonisation1 <- with(parms.temp, (del*rho1+(1-del)*Q.plus1)*(b-c1*rho1-c21*rho2-g*Q.plus2*p)*delta)
      #recolonisation2 <- with(parms.temp, (del*rho2+(1-del)*Q.plus2)*(b-c2*rho2-c12*rho1-g*(1-Q.plus1*n))*delta)
      # calculate recolonisation rates of all cells
      recolonisation1 <- with(parms.temp, (del*rho.nurse+(1-del)*Q.plus1)*(b-c1*Q.plus1-c21*Q.plus2-g*Q.plus2*p)*delta)
      recolonisation2 <- with(parms.temp, (del*rho.protege+(1-del)*Q.plus2)*(b-c2*Q.plus2-c12*Q.plus1-g*(1-Q.plus1*n))*delta)
      # calculate death rates                              
      death1 <- with(parms.temp, m*delta)
      death2 <- with(parms.temp, m*delta)
      
      # calculate regeneration rate and degradation rate
      regeneration <- with(parms.temp, (r + f*(Q.plus1+Q.plus2)*delta))
      degradation <- with(parms.temp, (d *delta))
      
      # check for sum of probabilities to be inferior 1 and superior 0
      if(any(c(recolonisation1+degradation,
               recolonisation2+degradation,  recolonisation1 + recolonisation2 + degradation,
               death1, death2, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in run",
                                                                  iteration, "time step", i, "! decrease delta!!!")) 
      if(any(c(recolonisation1, recolonisation2,
               degradation, regeneration,
               death1, death2) < 0)) warning(paste("a set probability falls below 0 in run",
                                                   iteration, "in time step", i, "! balance parameters!!!")) 
      # }}}
      # apply rules # {{{
      
      x.new$cells[which(x.old$cells == "0" & rnum <= recolonisation1)] <- "+1"
      x.new$cells[which(x.old$cells == "0" & rnum > recolonisation1 & rnum <= recolonisation1+recolonisation2)] <- "+2"
      x.new$cells[which(x.old$cells == "0" & rnum > recolonisation1 + recolonisation2 & rnum <= recolonisation1 + recolonisation2 + degradation)] <- "-"
      x.new$cells[which(x.old$cells == "+1" & rnum <= death1)] <- "0"
      x.new$cells[which(x.old$cells == "+2" & rnum <= death2)] <- "0"
      x.new$cells[which(x.old$cells == "-" & rnum <= regeneration)] <- "0"
      
      # 5 saving state of the new grid		
      
       
      result$rho.plus[i] <- sum( sum(x.new$cells == "+1"), sum(x.new$cells == "+2") )/(width*height)
      
      result$rho.nurse[i] <- sum(x.new$cells == "+1")/(width*height)
      
      result$rho.protege[i] <- sum(x.new$cells == "+2")/(width*height)
      
      result$q.1.1[i]  <- mean(subset(count(x.new, "+1")/4, x.new$cells == "+1") )

      result$q.2.2[i]  <- mean(subset(count(x.new, "+2")/4, x.new$cells == "+2") )
      
      result$q.1.2[i]  <- mean(subset(count(x.new, "+1")/4, x.new$cells == "+2") )

      result$q.2.1[i]  <- mean(subset(count(x.new, "+2")/4, x.new$cells == "+1") )
      
      
      x.old <- x.new
      
      
      if(i > t.min/delta+1 ) { #| result$rho.plus[i] == 0
     
        
        if( result$rho.plus[i] > 0) {
          t.1 <- (i-2*t.eval/delta):(i-t.eval/delta)-1
          t.2 <- (i-t.eval/delta):(i)
          stability.total <- (abs(mean(result$rho.plus[t.1]) - mean(result$rho.plus[t.2])))/(mean(result$rho.plus[t.1]))
          
          if( result$rho.nurse[i] > 0){
          stability.nurse <- (abs(mean(result$rho.nurse[t.1]) - mean(result$rho.nurse[t.2])))/(mean(result$rho.nurse[t.1]))
          } else {stability.nurse = 0}
          
          if( result$rho.protege[i] > 0){
          stability.protege <- (abs(mean(result$rho.protege[t.1]) - mean(result$rho.protege[t.2])))/(mean(result$rho.protege[t.1]))
          } else {stability.protege = 0}
        
          } else {
          stability.total = 0
          stability.nurse = 0
          stability.protege = 0
        }
        
        result$time[i] <- i*delta          
        
      }

    } # end of simulation run (over i)
    
    j = (j + 1)
    
    if(stability.total <= 0.0001 & stability.nurse <= 0.05 & stability.protege <= 0.05) { 
      nrep = (nrep + 1) 
      
      collect$lattice[[nrep]] <- x.new
      
      t.fin = result$time[length(result$time)]
      i.stable <- (length(result$time)-(2*t.eval/delta)):length(result$time)
      
      # switch for output at extinction
      if(result$rho.plus[t.fin] > 0) {
        result$out <- data.frame(
          
          replicate = j,
          
          rho.plus.ini = result$rho.plus[1],
          rho.plus = mean(result$rho.plus[i.stable], na.rm=T), 
          rho.plus.sd = sd(result$rho.plus[i.stable], na.rm=T),
          rho.plus.fin = result$rho.plus[length(result$time)],
          
          rho.nurse.ini = result$rho.nurse[1],
          rho.nurse = mean(result$rho.nurse[i.stable], na.rm=T), 
          rho.nurse.sd = sd(result$rho.nurse[i.stable], na.rm=T),
          rho.nurse.fin = result$rho.nurse[length(result$time)],
          
          rho.protege.ini = result$rho.protege[1],
          rho.protege = mean(result$rho.protege[i.stable], na.rm=T), 
          rho.protege.sd = sd(result$rho.protege[i.stable], na.rm=T),
          rho.protege.fin = result$rho.protege[length(result$time)],
          
          q.1.1.ini = result$q.1.1[1],
          q.1.1 = mean(result$q.1.1[i.stable], na.rm=T),
          q.1.1.sd = sd(result$q.1.1[i.stable], na.rm=T),
          q.1.1.fin = result$q.1.1[length(result$time)],
          
          q.2.2.ini = result$q.2.2[1],
          q.2.2 = mean(result$q.2.2[i.stable], na.rm=T),
          q.2.2.sd = sd(result$q.2.2[i.stable], na.rm=T),
          q.2.2.fin = result$q.2.2[length(result$time)],
          
          q.1.2.ini = result$q.1.2[1],
          q.1.2 = mean(result$q.1.2[i.stable], na.rm=T),
          q.1.2.sd = sd(result$q.1.2[i.stable], na.rm=T),
          q.1.2.fin = result$q.1.2[length(result$time)],
          
          q.2.1.ini = result$q.2.1[1],
          q.2.1 = mean(result$q.2.1[i.stable], na.rm=T),
          q.2.1.sd = sd(result$q.2.1[i.stable], na.rm=T),
          q.2.1.fin = result$q.2.1[length(result$time)],
          
          clus.1.1 = mean(result$q.1.1[i.stable]/result$rho.nurse[i.stable], na.rm=T),
          clus.1.1.sd = sd(result$q.1.1[i.stable]/result$rho.nurse[i.stable], na.rm=T),
          
          clus.2.2 = mean(result$q.2.2[i.stable]/result$rho.protege[i.stable], na.rm=T),
          clus.2.2.sd = sd(result$q.2.2[i.stable]/result$rho.protege[i.stable], na.rm=T),
          
          clus.1.2 = mean(result$q.1.2[i.stable]/result$rho.nurse[i.stable], na.rm=T),
          clus.1.2.sd = sd(result$q.1.2[i.stable]/result$rho.nurse[i.stable], na.rm=T),
          
          clus.2.1 = mean(result$q.2.1[i.stable]/result$rho.protege[i.stable], na.rm=T),
          clus.2.1.sd = sd(result$q.2.1[i.stable]/result$rho.protege[i.stable], na.rm=T),
          
          stability.total = stability.total,
          stability.nurse = stability.nurse,
          stability.protege = stability.protege,
          runtime = t.fin,
          starting = init.plant
        )
        
      } else {
        
        result$out <- data.frame(
          
          replicate = j,        
          
          rho.plus.ini = result$rho.plus[1],
          rho.plus = 0, 
          rho.plus.sd = 0,
          rho.plus.fin = 0,
          
          rho.nurse.ini = result$rho.nurse[1],
          rho.nurse = 0, 
          rho.nurse.sd = 0,
          rho.nurse.fin = 0,
          
          rho.protege.ini = result$rho.protege[1],
          rho.protege = 0, 
          rho.protege.sd = 0,
          rho.protege.fin = 0,
          
          q.1.1.ini = result$q.1.1[1],
          q.1.1 = NA,
          q.1.1.sd = NA,
          q.1.1.fin = NA,
          
          q.2.2.ini = result$q.2.2[1],
          q.2.2 = NA,
          q.2.2.sd = NA,
          q.2.2.fin = NA,
          
          q.1.2.ini = result$q.1.2[1],
          q.1.2 = NA,
          q.1.2.sd = NA,
          q.1.2.fin = NA,
          
          q.2.1.ini = result$q.2.1[1],
          q.2.1 = NA,
          q.2.1.sd = NA,
          q.2.1.fin = NA,
          
          clus.1.1 = NA,
          clus.1.1.sd = NA,
          
          clus.2.2 = NA,
          clus.2.2.sd = NA,
          
          clus.1.2 = NA,
          clus.1.2.sd = NA,
          
          clus.2.1 = NA,
          clus.2.1.sd = NA,
          
          
          stability.total = 0,
          stability.nurse = 0,
          stability.protege = 0,
          runtime = t.fin,
          starting = init.plant
        ) 
      }
      
      collect$out[nrep,] <- result$out
    }
   
    
  }#-> collect$out
  
  
  result <- list()
  
  # pool replicates to one row (means, sd)
  stable <- collect$out$stability.total < 0.0001 & collect$out$stability.nurse < 0.05 & collect$out$stability.protege < 0.05
  
  result$runs <- collect$out 
  result$grids <- collect$lattice
  
  
  result$out <- data.frame(
    ID = iteration,
    del = parms.temp$del,
    g = parms.temp$g, 
    b = parms.temp$b, 
    m = parms.temp$m,
    com = parms.temp$com,
    starting = parms.temp$starting,
    c1 = parms.temp$c1,
    c2 = parms.temp$c2,
    c12 = parms.temp$c12,
    c21 = parms.temp$c21,
    n = parms.temp$n, 
    f = parms.temp$f, 
    r = parms.temp$r, 
    d = parms.temp$d, 
    d = parms.temp$d, 
    
    nrep = length(which(stable)),
    runtime = mean(result$runs$runtime[stable], na.rm=T), 
    runtime.sd = sd(result$runs$runtime[stable], na.rm=T), 
    
    rho.plus.ini = mean(result$runs$rho.plus.ini[stable], na.rm=T),
    rho.plus = mean(result$runs$rho.plus[stable], na.rm=T),
    rho.plus.sd = mean(result$runs$rho.plus.sd[stable], na.rm=T),
    rho.plus.sd.inter.replicates = sd(result$runs$rho.plus[stable], na.rm=T),
    
    rho.nurse.ini = mean(result$runs$rho.nurse.ini[stable], na.rm=T),
    rho.nurse = mean(result$runs$rho.nurse[stable], na.rm=T), 
    rho.nurse.sd = mean(result$runs$rho.nurse.ini[stable], na.rm=T),
    rho.nurse.sd.inter.replicates = sd(result$runs$rho.nurse[stable], na.rm=T),
    
    rho.protege.ini = mean(result$runs$rho.protege.ini[stable], na.rm=T),
    rho.protege = mean(result$runs$rho.protege[stable], na.rm=T), 
    rho.protege.sd = mean(result$runs$rho.protege.ini[stable], na.rm=T),
    rho.nurse.sd.inter.replicates = sd(result$runs$rho.protege[stable], na.rm=T),
    
    q.1.1 =  mean(result$runs$q.1.1[stable], na.rm=T),
    q.1.1.sd = mean(result$runs$q.1.1.sd[stable], na.rm=T),
    q.1.1.inter.replicates = sd(result$runs$q.1.1[stable], na.rm=T),
    
    q.2.2 =  mean(result$runs$q.2.2[stable], na.rm=T),
    q.2.2.sd = mean(result$runs$q.2.2.sd[stable], na.rm=T),
    q.2.2.inter.replicates = sd(result$runs$q.2.2[stable], na.rm=T),
    
    q.1.2 =  mean(result$runs$q.1.2[stable], na.rm=T),
    q.1.2.sd = mean(result$runs$q.1.2.sd[stable], na.rm=T),
    q.1.2.inter.replicates = sd(result$runs$q.1.2[stable], na.rm=T),
    
    q.2.1 =  mean(result$runs$q.2.1[stable], na.rm=T),
    q.2.1.sd = mean(result$runs$q.2.1.sd[stable], na.rm=T),
    q.2.1.inter.replicates = sd(result$runs$q.2.1[stable], na.rm=T),
    
    clus.1.1 = mean(result$runs$clus.1.1[stable], na.rm=T),
    clus.1.1.sd = sd(result$runs$clus.1.1.sd[stable], na.rm=T),
    clus.1.1.inter.replicates = sd(result$runs$clus.1.1[stable], na.rm=T),
    
    clus.2.2 = mean(result$runs$clus.2.2[stable], na.rm=T),
    clus.2.2.sd = sd(result$runs$clus.2.2.sd[stable], na.rm=T),
    clus.2.2.inter.replicates = sd(result$runs$clus.2.2[stable], na.rm=T),
    
    clus.1.2 = mean(result$runs$clus.1.2[stable], na.rm=T),
    clus.1.2.sd = sd(result$runs$clus.1.2.sd[stable], na.rm=T),
    clus.1.2.inter.replicates = sd(result$runs$clus.1.2[stable], na.rm=T),
    
    clus.2.1 = mean(result$runs$clus.2.1[stable], na.rm=T),
    clus.2.1.sd = sd(result$runs$clus.2.1.sd[stable], na.rm=T),
    clus.2.1.inter.replicates = sd(result$runs$clus.2.1[stable], na.rm=T),
    
    runtime = t.fin
  )
  
  
  
  
  save(result, file = paste("result.coex.17.06", iteration,".rdata", sep = ""))
  
  gc() 
  return(result$out)
  
}-> output

write.csv(output, file = paste("output.coex.17.06.csv", sep = "") )# }}}

stopCluster(cl)
