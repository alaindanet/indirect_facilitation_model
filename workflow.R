library(simecol)
library(tidyverse)
library(magrittr)

devtools::document()
#devtools::use_vignette("three_states_model")

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

parms(mod)["g"] <- 0.13
parms(mod)["gamma1"] <- 0.10
parms(mod)["protection_type"] <- list("first_protect")
parms(mod)["protection_type"] <- list("linear")
times(mod) <- c(from = 0, to = 3000, by = 1)
mod_run <- sim(mod)
tail(out(mod_run), 40)
plotnp(mod_run)

run_2d_gradient

################################
#  Grazing and gamma gradient  #
################################

# Run_gradient g and gamma1
g_gradient <- seq(0, 0.3, length.out = 10)
gamma_gradient <- seq(0, 0.3, length.out = 10)
gradient_2d <- run_2d_gradient(
  y = "g",
  gradienty = g_gradient,
  x = "gamma1",
  gradientx = gamma_gradient,
  param = c(b = 0.8, u = 0, tau_n = 0, protection_type = list("linear")),
  time_seq = c(from = 0, to = 3000, by = 1),
  nb_cores = 4,
  solver_type = NULL#steady_state
  )

gradient_2d
averaged_runs <- avg_runs(gradient_2d, cut_row = 10)

plot_diagram(averaged_runs, param = c(x = "gamma1", y = "g"), debug_mode = FALSE)
#ggsave("diag_g_gamma_linear_tau_n=0.png")
plot_diagram(averaged_runs, debug_mode = TRUE)
# Not work bc many values by gradient value 

##################################
#  Grazing and aridity gradient  #
##################################

g_gradient <- seq(0, 0.3, length.out = 10)
b_gradient <- seq(1, 0, length.out = 10)
gradient_2d <- run_2d_gradient(
  y = "g",
  gradienty = g_gradient,
  x = "b",
  gradientx = b_gradient,
  param = c(protection_type = list("first_protect"), gamma1 = 0.1, u = 5),
  time_seq = c(from = 0, to = 5000, by = 1),
  nb_cores = 4,
  solver_type = NULL#steady_state
  )

averaged_runs <- avg_runs(gradient_2d, cut_row = 10)

plot_diagram(averaged_runs, param = c(x = "b", y = "g"))
ggsave("diag_aridity_grazing_first_protect_u=5.png")
