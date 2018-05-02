
upca <- upca_model()
parms(upca)["wstar"] <- 0.1
equations(upca)$f <- equations(upca)$f1
solver(upca) <- steady_state_upca
test <- sim(upca)
plotupca(test)

# Test
mod <- indirect_facilitation_model()
times(mod) <- c(from = 0, to = 1000, by = .3)
solver(mod) <- "lsoda"
mod_run <- sim(mod)
