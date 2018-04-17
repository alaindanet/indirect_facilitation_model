library(simecol)
library('tidyverse')
library(magrittr)

devtools::document()
#devtools::use_vignette("three_states_model")

upca <- upca_model()
parms(upca)["wstar"] <- 0.1
equations(upca)$f <- equations(upca)$f1
solver(upca) <- steady_state_upca
test <- sim(upca)
plotupca(test)

root <- function(time, init, parms) {
  dstate <- unlist(upca_ode(time, init, parms))
  return(sum(abs(dstate)) - 1e-4)
}
lsodar(init(upca), c(0, 1e10), upca_ode, parms(upca), rootfun = root)

##   

# Test
mod <- indirect_facilitation_model()
mod
times(mod) <- c(from = 0, to = 1000, by = .3)
solver(mod) <- "lsoda"
mod_run <- sim(mod)

parms(mod)["g"] <- 0.03
parms(mod)["gamma1"] <- 0.10
times(mod) <- c(from = 0, to = 1000, by = 1)
mod_run <- sim(mod)
#tail(out(mod_run), 40)
plotnp(mod_run)

# Run_gradient
g_gradient <- seq(0, 0.2, 0.02)
gamma_gradient <- seq(0, 0.2, length.out = 10)

spec_param <- parms(mod)["b"] <- .8

b_gradient <- seq(1, 0, length.out = 10)
gradient_2d <- run_2d_gradient(
  y = "g",
  gradienty = g_gradient,
  x = "b",
  gradientx = b_gradient,
  param = c(gamma1 = 0.1),
  time_seq = c(from = 0, to = 5000, by = 1),
  nb_cores = 4,
  solver_type = NULL#steady_state
  )

gradient_2d
averaged_runs <- gradient_2d %>%
  dplyr::mutate(avg = purrr::map(runs, avg_runs, cut_row = 10)) %>%
  tidyr::unnest(avg)

plot_diagram(averaged_runs, param = c(x = "b", y = "g"))
#ggsave("diag_aridity_grazing.png")
plot_diagram(averaged_runs, debug_mode = TRUE)
# Not work bc many values by gradient value 
plotnp_gradient(averaged_runs)

# filter
test <- dplyr::filter(averaged_runs, g == 0.02 & gamma1 == 0.02)
test2 <- dplyr::filter(gradient_2d, g == 0 & gamma1 == 0)
tail(test[1, ]$runs[[1]])
plot(test2[1, ]$runs[[1]]$N)

