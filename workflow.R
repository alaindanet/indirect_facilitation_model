library(simecol)
library(tidyverse)
library(magrittr)

devtools::document()
#devtools::use_vignette("three_states_model")

mod <- indirect_facilitation_model()
init(mod)
parms(mod)["g"] <- 0.13
parms(mod)["gamma1"] <- 0.1
parms(mod)["protection_type"] <- list("first_protect")
parms(mod)["protection_type"] <- list("linear")
times(mod) <- c(from = 0, to = 3000, by = 1)
mod_run <- sim(mod)
tail(out(mod_run), 40)
plotnp(mod_run)

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

###############################
#  Bifurcation state diagram  #
###############################

bifurc <- run_bifurcation(
  gradientx = seq(0.2, 1, by = 0.1),
  gradienty = c(0.4, 0.01),
  param = c(
    g = .20, gamma1 = .15,
    protection_type = list("first_protect"), #list("first_protect")
    u = 5
    )
  )

averaged_runs <- avg_runs(bifurc, cut_row = 1)
plotnp(averaged_runs, alpha = 0.65)


################################################################################
#                              Four states model                               #
################################################################################

mod <- two_facilitation_model()
init(mod)
parms(mod)["g"] <- 0.05
parms(mod)["gamma1"] <- 0.01
parms(mod)["protection_type"] <- list("first_protect")
parms(mod)["protection_type"] <- list("linear")
times(mod) <- c(from = 0, to = 3000, by = 1)
mod_run <- sim(mod)
tail(out(mod_run), 40)
plotnp(mod_run)

##################################
#  Grazing and aridity gradient  #
##################################

g_gradient <- seq(0, 0.3, length.out = 10)
b_gradient <- seq(1, .3, length.out = 10)
gradient_2d <- run_2d_gradient(
  y = "g",
  gradienty = g_gradient,
  x = "b",
  gradientx = b_gradient,
  model_spec = two_facilitation_model(),
  param = c(protection_type = list("linear"), gamma1 = 0.1, n = 1),
  time_seq = c(from = 0, to = 3000, by = 1),
  nb_cores = 4,
  solver_type = NULL#steady_state
  )

averaged_runs <- avg_runs(gradient_2d, cut_row = 10)
plot_diagram(averaged_runs, param = c(x = "b", y = "g"))
#ggsave("diag_aridity_grazing_first_protect_u=5.png")


