#'#'#############################################################
#'  Analyse simulation results of pair-approximations models
#'  Author: Alain Danet
#'  
#'#############################################################'

setwd("~/git/thesis_chap1_model/simulation/model2/")

data <- read.table("output.csv", sep=",", header=T, row.names=1)

# plot relationship between species co-occurences and others parameters

pdf("Model2_CnpVSparameters.pdf")
par(mfrow=c(3,3))
plot(Cnp~g + b + cn + cp + cpn + cnp + n + del,data)
dev.off()
par(mfrow=c(1,1))

# Simulations with positive co-occurences
gg <- data[data$Cnp >=1.0 & data$del==0.9,]
par(mfrow=c(3,3))
plot(Cnp~g + b + cn + cp + cpn + cnp + n + del,gg)
par(mfrow=c(1,1))
dev.off()

par(mfrow=c(2,1))
plot(Cnp~rhon, gg, col= rainbow(20)[ as.factor(round(gg$rhop,1)) ] )
legend(x=0.1,
       y=2.0,
       legend= levels(as.factor(round(gg$rhop,1))),
       col = rainbow(20),
       pch=1)
plot(Cnp~rhop, gg, col= rainbow(20)[ as.factor(round(gg$rhon,1)) ] )
legend(x=0.1,
       y=2.0,
       legend= levels(as.factor(round(gg$rhon,1))),
       col = rainbow(20),
       pch=1)
