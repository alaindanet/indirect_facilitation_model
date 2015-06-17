#'#'#####################################
#'  Function to analyses simulations
#'  Author: Alain Danet
#'#####################################'

# Pair simulation# {{{

# Function to plot coexistence diagrams pair simulation# {{{
diag.pairs <- function(data  = data,
		      delta = 0.1,
		      cp    = 0.2,
		      cn    = 0.2,
		      cpn   = 0.15,
		      cnp   = 0.1,
		      m     = 0.10,
		      n     = 1,
		      alpha = 150,
		      main = "Coexistence"){

	data <- data[data$del==delta & data$cp==cp & data$cn==cn & round(data$cpn,2)==cpn & round(data$cnp,2)==cnp & round(data$m,2)==m & data$n==n,]

	b=unique(data$b)
	g=unique(data$g)

	if(nrow(data)==0) warning (" Non valid parameters combination")
	if(!is.vector(b) | !is.vector(g)) warning ("b & g must be a numeric vector of unique values of each parameter")
	mat1 <- matrix(NA,
		       ncol=length(g),
		       nrow=length(b)
		       )
	mat2 <- mat1
	mat3 <- mat1

	for(i in 1:nrow(data)){
		mat1[which(unique(b)== data[i,]$b), which(unique(g)== data[i,]$g)] <- ifelse(data[i,]$rhon >= 0.05, 1, NA)
		mat2[which(unique(b)== data[i,]$b), which(unique(g)== data[i,]$g)] <- ifelse(data[i,]$rhop >= 0.05, 1, NA)
		mat3[which(unique(b)== data[i,]$b), which(unique(g)== data[i,]$g)] <- ifelse(data[i,]$rhop <= 0.05 & data[i,]$rhon <= 0.05 , 1, NA)
	}

	green <- col2rgb(c("green1")) # Sp1
	yellow <- col2rgb(c("yellow2")) # Desert
	grey <- col2rgb("black") # Sp2
 	b  <- 1-b

	image(b, g, mat2, col= rgb( grey[1,], grey[2,], grey[3,], maxColorValue = 255, alpha = 255), xlab="Aridité", ylab="Pâturage", main=main)
	image(b, g, mat1, col= rgb( green[1,], green[2,], green[3,], maxColorValue = 255, alpha = alpha), add=T)
# 	image(b, g, mat3, col= rgb( yellow[1,], yellow[2,], yellow[3,], maxColorValue = 255, alpha = alpha), add=T)
}# }}}

# Function to plot co-occurence diagrams # {{{
diag.c.pairs <- function(data  = data,
		      delta = 0.1,
		      cp    = 0.2,
		      cn    = 0.2,
		      cpn   = 0.15,
		      cnp   = 0.1,
		      m     = 0.10,
		      n     = 1,
		      alpha = 150,
		      main = "Co-occurences"){

	data <- data[data$del==delta & data$cp==cp & data$cn==cn & round(data$cpn,2)==cpn & round(data$cnp,2)==cnp & round(data$m,2)==m & data$n==n,]

	b=unique(data$b)
	g=unique(data$g)

	if(nrow(data)==0) warning (" Non valid parameters combination")
	if(!is.vector(b) | !is.vector(g)) warning ("b & g must be a numeric vector of unique values of each parameter")
	mat <- matrix(NA,
		       ncol=length(g),
		       nrow=length(b)
		       )

	for(i in 1:nrow(data)){
		mat[which(unique(b)== data[i,]$b), which(unique(g)== data[i,]$g)] <- data[i,]$Cnp
	}

 	b  <- 1-b

	image(b, g, mat, xlab="Aridité", ylab="Pâturage", main=main, zlim=c(0,1.2))
}# }}}

# Function to plot the dynamics of the system along an aridity gradient.# {{{
stab.pairs <- function(data  = data,
		       delta = 0.1,
		       cp    = 0.2,
 		       cn    = 0.2,
 		       cpn   = 0.15,
 		       cnp   = 0.1,
 		       m     = 0.10,
 		       n     = 1,
 		       g     = 0.0,
 		       alpha = 150,
 		       main = "Dynamic"){

	data <- data[data$del==delta & data$cp==cp & data$cn==cn & round(data$cpn,2)==cpn & round(data$cnp,2)==cnp & round(data$m,2)==m & data$n==n & round(data$g,2)==g,]

	if(nrow(data)==0) warning (" Non valid parameters combination")

	data$b = 1-data$b
	plot(rhop~b,
	     data,
	     type="l",
	     col="black",
	     xlab="Aridité",
	     ylab="Rho et Clustering",
	     ylim=c(0,1.5),
	     main=main)
	lines(rhon~b,data, col="green")
	lines(Cpp~b,data,lty=2, col="black")
	lines(Cnn~b,data,lty=2, col="green")
	lines(Cnp~b,data,lty=2, col="red")

}# }}}# }}}

# Cellular automata simulation# {{{

diag.coex <- function(folder          = "~/result/scenario1.1/",
		      name            = "result_coex_scenario",
		      file_iterations = "CA_output.csv",
		      com             = "poly",
		      b               = seq(0.7,0.5,-0.1),
		      g               = seq(0,0.1,0.05),
                      delta           = 0.1,
		      n               = 1,
		      c_12            = 0.05,
		      c_21            = 0.05,
		      m               = 0.2,
		      c_1             = 0.2,
		      c_2             = 0.2,
		      seuil           = 0.9,
		      alpha           = 150,
		      main            = "Some title", ...){
  
  # Read iteration table
  iterations <- read.table(paste(folder, file_iterations, sep=""), sep=",", dec=".", row.names=1, header=T)
  
  # Load all iteration ID for choosen variable combination
  iter <- iterations[iterations$com==com  & iterations$n == n & iterations$c_12 == c_12 & iterations$c_21 == c_21, ]$ID 
  if(length(iter)==0) warning (" Non valid parameters combination")
  if(!is.vector(iter)) warning ("iter must be a numeric vector of iterations ID")
  if(!is.vector(b) | !is.vector(g)) warning ("b & g must be a numeric vector of unique values of each parameter")
  
  # set up matrices to fill result: SP1 wins, SP2 wins or no one  is alive along aridity (b) and grazing (g) gradients.
  mat1 <- matrix(NA,ncol=length(g),nrow=length(b))
  mat2 <- mat1
  mat3 <- mat1
  for(n in 1:length(iter)){
    infile <- paste(folder,name,iter[n],".rdata",sep="")# load each simulation of choosen iteration
    load(infile)
   
    parms_temp <- result$out[,c("g","b","clus_1_2","clus_1_1", "clus_2_2","clus_2_1")] # take parameters
    parms_temp$g  <- round(parms_temp$g,2)
    # Fill the matrices
    mat1[which(unique(b)==parms_temp$b), which(unique(g)==parms_temp$g) ] <- ifelse( mean(result$runs$rho_nurse > 0.05) >= seuil, 1, NA) 
    mat2[which(unique(b)==parms_temp$b), which(unique(g)==parms_temp$g) ] <- ifelse( mean(result$runs$rho_protege > 0.05) >= seuil, 2, NA)
    mat3[which(unique(b)==parms_temp$b), which(unique(g)==parms_temp$g)] <- ifelse( mean(result$runs$rho_protege < 0.05) >= seuil & mean(result$runs$rho_nurse < 0.05) >= seuil, 3, NA)
  }
  
  green <- col2rgb(c("green1")) # Sp1
  yellow <- col2rgb(c("yellow2")) # Desert
  grey <- col2rgb("black") # Sp2
  
  b <- 1 - unique(iterations$b)
  
  # Plot the matrices 
  if(com == "poly"){
    image(b, g ,mat2,col=rgb( grey[1,], grey[2,], grey[3,], maxColorValue = 255, alpha = 255),
          xlab="Aridité", ylab="Pâturage", main=main)
    image(b , g, mat1, col = rgb( green[1,], green[2,], green[3,], maxColorValue = 255, alpha = alpha), add=T)
    #image(b,g,mat3,col=rgb( blue[1,], blue[2,], blue[3,], maxColorValue = 255, alpha = alpha), add=T)
    image(b , g, mat3, col = rgb( yellow[1,], yellow[2,], yellow[3,], maxColorValue = 255, alpha = alpha), add=T)
  }
  if(com == "mono_p"){
    image(b,g,mat2,col=rgb( grey[1,], grey[2,], grey[3,], maxColorValue = 255, alpha = alpha), family="Times", ...)
  }
  if(com == "mono_n"){
    image(b , g, mat1, col = rgb( green[1,], green[2,], green[3,], maxColorValue = 255, alpha = alpha), family="Times", ...)
  }

}# }}}

# Cellular automata 2: replacement of _ by . #{{{

diag.coex2 <- function(folder          = "~/result/scenario1.1/",
		      name            = "result.coex.scenario",
		      file.iterations = "CA.output.csv",
		      com             = "poly",
		      b               = seq(0.7,0.5,-0.1),
		      g               = seq(0,0.1,0.05),
                      delta           = 0.1,
		      n               = 1,
		      c12            = 0.05,
		      c21            = 0.05,
		      m               = 0.2,
		      c1             = 0.2,
		      c2             = 0.2,
		      seuil           = 0.9,
		      alpha           = 150,
		      main            = "Some title", ...){
  
  # Read iteration table
  iterations <- read.table(paste(folder, file.iterations, sep=""), sep=",", dec=".", row.names=1, header=T)
  
  # Load all iteration ID for choosen variable combination
  iter <- iterations[iterations$com==com  & iterations$n == n & iterations$c12 == c12 & iterations$c21 == c21, ]$ID 
  if(length(iter)==0) warning (" Non valid parameters combination")
  if(!is.vector(iter)) warning ("iter must be a numeric vector of iterations ID")
  if(!is.vector(b) | !is.vector(g)) warning ("b & g must be a numeric vector of unique values of each parameter")
  
  print(length(iter))
  # set up matrices to fill result: SP1 wins, SP2 wins or no one  is alive along aridity (b) and grazing (g) gradients.
  mat1 <- matrix(NA,ncol=length(g),nrow=length(b))
  mat2 <- mat1
  mat3 <- mat1
  for(n in 1:length(iter)){
    infile <- paste(folder,name,iter[n],".rdata",sep="")# load each simulation of choosen iteration
    load(infile)
   
    parms.temp <- result$out[,c("g","b","clus.1.2","clus.1.1", "clus.2.2","clus.2.1")] # take parameters
    parms.temp$g  <- round(parms.temp$g,2)
    # Fill the matrices
    mat1[which(unique(b)==parms.temp$b), which(unique(g)==parms.temp$g) ] <- ifelse( mean(result$runs$rho.nurse > 0.05) >= seuil, 1, NA) 
    mat2[which(unique(b)==parms.temp$b), which(unique(g)==parms.temp$g) ] <- ifelse( mean(result$runs$rho.protege > 0.05) >= seuil, 2, NA)
    mat3[which(unique(b)==parms.temp$b), which(unique(g)==parms.temp$g)] <- ifelse( mean(result$runs$rho.protege < 0.05) >= seuil & mean(result$runs$rho.nurse < 0.05) >= seuil, 3, NA)
  }
  print(mat2)
  green <- col2rgb(c("green1")) # Sp1
  yellow <- col2rgb(c("yellow2")) # Desert
  grey <- col2rgb("black") # Sp2
  
  b <- 1 - unique(iterations$b)
  
  # Plot the matrices 
  if(com == "poly"){
    image(b, g ,mat2,col=rgb( grey[1,], grey[2,], grey[3,], maxColorValue = 255, alpha = 255),
          xlab="Aridité", ylab="Pâturage", main=main)
    image(b , g, mat1, col = rgb( green[1,], green[2,], green[3,], maxColorValue = 255, alpha = alpha), add=T)
    #image(b,g,mat3,col=rgb( blue[1,], blue[2,], blue[3,], maxColorValue = 255, alpha = alpha), add=T)
    image(b , g, mat3, col = rgb( yellow[1,], yellow[2,], yellow[3,], maxColorValue = 255, alpha = alpha), add=T)
  }
  if(com == "mono.p"){
    image(b,g,mat2,col=rgb( grey[1,], grey[2,], grey[3,], maxColorValue = 255, alpha = alpha), family="Times", ...)
  }
  if(com == "mono.n"){
    image(b , g, mat1, col = rgb( green[1,], green[2,], green[3,], maxColorValue = 255, alpha = alpha), family="Times", ...)
  }

}
# }}}
