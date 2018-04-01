
devtools::document()
library(simecol)

upca <- upca_model()
parms(upca)
equations(upca)$f <- equations(upca)$f2
test <- sim(upca)
plotupca(test)

##   

# Test
mod <- indirect_facilitation_model()
mod
parms(mod)
times(mod) <- c(from = 0, to = 1000, by = 10)
equations(mod)$asymp_cost <- equations(mod)$asymp_cost
test <- sim(mod)
out(test)
plot(test)
traceback()
