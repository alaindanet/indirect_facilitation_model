library(simecol)

devtools::document()

upca <- upca_model()
parms(upca)
equations(upca)$f <- equations(upca)$f2
test <- sim(upca)
plotupca(test)

##   

# Test
mod <- indirect_facilitation_model()
mod
times(mod) <- c(from = 0, to = 10000, by = 100)
equations(mod)$asymp_cost <- equations(mod)$asymp_cost
mod_run <- sim(mod)

parms(mod)["g"] <- 0.06291
mod_run <- sim(mod)
plotnp(mod_run)



g_gradient <- seq(0, 0.2, length.out = 20)
gamma_gradient <- seq(0, 0.2, length.out = 20)

test  <- run_2d_gradient()


library(tibble)
library(magrittr)
test
test[1,]$runs[[1]][1,]

averaged <- test %>%
  dplyr::mutate(avg = purrr::map(runs, avg_runs, cut_row = 100)) %>%
  tidyr::unnest(avg)
