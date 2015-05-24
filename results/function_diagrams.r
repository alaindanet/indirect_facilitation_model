#################################################################################
#
# Function to plot coexistence diagrams
# Author: Alain Danet
#
# Date: 7.01.2015
#
#################################################################################


# Tasks to do:
# ============
# Thinking to a smart way to define a general function

# It plots a diagram of survival species across aridity and grazing gradient

# Desctiption of the variables:
# ============================
# folder: path to the .rdata objects
# name: common name part of .rdata objects
# iterations: output .csv with summary results and interation IDs
# com: can be poly, mono_p, mono_n. Depends if simulations are mono or polyspecific.
# ...Others parameters of the model
# seuil: threshold proportion of replicates which indicates same result to take in account this result
# alpha: transparency factor of colours

diag.coex4 <- function(folder="~/result/scenario1.1/", name="result_scenario1.1_", file_iterations="output_scenario1.1.csv",com= "poly", b=seq(0.7,0.5,-0.1), g=seq(0,0.1,0.05),
                       delta=0.1, n = 1, c_12 = 0.05,c_21 = 0.05, m=0.2, c_1 = 0.2, c_2 = 0.2, seuil=0.9,  alpha=150, ...){
  
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
          xlab="Aridité", ylab="Pâturage")#, family="Times")
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

}

