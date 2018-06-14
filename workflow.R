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
  param = c(protection_type = list("first_protect"), gamma1 = 0.1, u = 0),
  time_seq = c(from = 0, to = 5000, by = 1),
  nb_cores = 4,
  solver_type = NULL#steady_state
  )

averaged_runs <- avg_runs(gradient_2d, cut_row = 10)
plot_diagram(averaged_runs, param = c(x = "b", y = "g"))
ggsave("./inst/figs/diag_aridity_grazing_first_protect_u=0.png")

###############################
#  Bifurcation state diagram  #
###############################

bifurc <- run_bifurcation(
  gradientx = seq(0.2, 1, by = 0.1),
  gradienty = c(0.4, 0.01),
  param = c(
    g = .15, gamma1 = .1,
    protection_type = list("first_protect"), #list("first_protect")
    u = 0 
    )
  )

averaged_runs <- avg_runs(bifurc, cut_row = 1)
plotnp(averaged_runs, alpha = 0.65)


################################################################################
#                              Four states model                               #
################################################################################

mod <- two_facilitation_model()
init(mod)
names(init(mod))
parms(mod)["g"] <- 0.3
parms(mod)["gamma1"] <- 0.08
parms(mod)["protection_type"] <- list("first_protect")
parms(mod)["protection_type"] <- list("linear")
times(mod) <- c(from = 0, to = 1000, by = 1)
parms(mod)["u"] <- 30
mod_run <- sim(mod)
plotnp(mod_run)

#############################################
#  Grazing and grazing protection strength  #
#############################################

g_gradient <- seq(0, 0.3, length.out = 10)
u_gradient <- seq(0, 20, length.out = 10)
gradient_2d <- run_2d_gradient(
  y = "g",
  gradienty = g_gradient,
  x = "u",
  gradientx = u_gradient,
  model_spec = "two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1),
  time_seq = c(from = 0, to = 3000, by = 1),
   nb_cores = 4,
  solver_type = NULL
  )

averaged_runs <- avg_runs(gradient_2d, cut_row = 10)
plot_diagram(averaged_runs, param = c(x = "u", y = "g"))
ggsave("inst/figs/four_states/diag_u_grazing_first_protect_gamma1=.1.png")

##################################
#  Grazing and aridity gradient  #
##################################

g_gradient <- seq(0, 0.3, length.out = 30)
b_gradient <- seq(1, 0.2, length.out = 30)
gradient_2d <- run_2d_gradient(
  y = "g",
  gradienty = g_gradient,
  x = "b",
  gradientx = b_gradient,
  model_spec = "two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1, u = 5),
  time_seq = c(from = 0, to = 3000, by = 1),
  nb_cores = 4,
  solver_type = NULL
  )

averaged_runs <- avg_runs(gradient_2d, cut_row = 10)
plot_diagram(averaged_runs, param = c(x = "b", y = "g"))
#save(gradient_2d, file = "diag_aridity_grazing_first_protect_u=5.RData" )
ggsave("inst/figs/four_states/diag_aridity_grazing_first_protect_u=5.pdf",
  scale = .8)

###############################
#  Bifurcation state diagram  #
###############################

bifurc <- run_bifurcation(
  gradientx = seq(0.2, 1, length.out = 30),
  gradienty = c(0.4, 0.01),
  model_spec = "two_facilitation_model",
  time_seq = c(from = 0, to = 3000, by = 1),
  param = c(
    g = .3, gamma1 = .1,
    protection_type = list("first_protect"), #list("first_protect")
    u = 5
    )
  )

averaged_runs <- avg_runs(bifurc, cut_row = 1)
plotnp(averaged_runs, alpha = 0.65)
ggsave("inst/figs/four_states/bifurc_first_protect_u=0_g=.25_gamma1=.1.png")

##############
#  Scenarii  #
##############

g_gradient <- seq(0, 0.3, length.out = 30)
b_gradient <- seq(1, 0.2, length.out = 30)

u_test <- c(0, 5)
output <- vector(mode = "list", length = 3)
names(output) <- sapply(u_test, function (x){paste("u", x, sep = "") })

for (i in seq_along(u_test)) {
  output[[i]] <- run_scenarii_gradient(
    y = "g",
    gradienty = g_gradient,
    x = "b",
    gradientx = b_gradient,
    model_spec = "two_facilitation_model",
    param = c(protection_type = list("first_protect"), gamma1 = 0.1, u = u_test[i]),
    time_seq = c(from = 0, to = 3000, by = 1),
    nb_cores = 3,
    solver_type = NULL,
    scenarii = init_scenarii(type = "bifurcation", ini_cover = .8, low_cover = .01)
    )
}
save(output, file = "./inst/scenar_bifurc_u=0_5_gamma1_.1.Rdata")

scenar_avg <- avg_runs(scenar)

g <- plot_diagram(scenar_avg, param = c(x = "b", y = "g", type = "scenario"),
  debug_mode = FALSE)
g

test <- compute_states(scenar_avg, param = c(x = "b", y = "g", type = "scenario"))
test
test$run %>%
  spread(scenario, state)


plotnp(scenar_avg)
ggsave("inst/figs/four_states/scenar_bifurc_first_protect_u=0_test.png")

names(output)

test <- avg_runs(output$u10)
plot_diagram(test)

##############
#  Figure 1  #
##############
library("stringr")

load(file = "./inst/scenar_bifurc_u=0_10_15_gamma1_.1.Rdata")
test <- lapply(output, FUN = avg_runs)
rm(output)
test2 <- lapply(test, FUN = compute_states, param = c(x = "b", y = "g", type = "scenario"))

test3 <- list()
for (i in seq_along(names(test2))) {
  test3[[i]] <- test2[[names(test2)[i]]] %>%
    dplyr::mutate(
      u = as.numeric(substr(names(test2)[i], 2, 10)),
      state = stringr::str_replace(state, "warning", "extinct")
      )
}
allu <- bind_rows(test3)

# À gauche 
u0 <- allu %>%
  filter(scenario == "together", u == 0, b >= 0.5)
class(u0) <- c("tibble", "data.frame", "states")
plot_diagram(u0, param = c(x = "b", y = "g"))
ggsave("inst/figs/four_states/diag_u=0.pdf", width = 7, height = 5, units =
  "cm")

# À droite 
load("diag_aridity_grazing_first_protect_u=5.RData") # gradient_2d

u5 <- avg_runs(gradient_2d)

compute_states(u5, param = c(x = "b", y = "g"))$run

states_u5 <- compute_states(u5, param = c(x = "b", y = "g")) %>%
  filter(b >= 0.5)
class(states_u5) <- c("tibble", "data.frame", "states")
plot_diagram(states_u5, param = c(x = "b", y = "g"))
ggsave("inst/figs/four_states/diag_u=5.pdf", width = 7, height = 5, units = "cm")


# Bonus
u10 <- allu %>%
  filter(scenario == "together", u == 10)
class(u10) <- c("tibble", "data.frame", "states")
plot_diagram(u10, param = c(x = "b", y = "g"))
ggsave("inst/figs/four_states/diag_u=10.pdf", scale = .8)

##############
#  Figure 2  #
##############

#u_grad <- c(0, 5)
u_grad <- c(0)
#g_grad <- c(0, .05, .1, .2, .3)
g_grad <- c(.09)

for(u in seq_along(u_grad)){
  for(g in seq_along(g_grad)) {
  
    bifurc <- run_bifurcation(
      gradientx = seq(0.2, 1, length.out = 30),
      gradienty = c(0.4, 0.01),
      model_spec = "two_facilitation_model",
      time_seq = c(from = 0, to = 3000, by = 1),
      param = c(
	g = g_grad[g], gamma1 = .1,
	protection_type = list("first_protect"),
	u = u_grad[u]
	)
      )

    averaged_runs <- avg_runs(bifurc, cut_row = 1)
    plotnp(averaged_runs, alpha = 0.65)
    plot_name <- paste("inst/figs/four_states/bifurc_first_protect_u=",
      u_grad[u], "_g=", g_grad[g], "_gamma1=.1.pdf", sep = "")
    ggsave(plot_name, width = 7, height = 5, units = "cm")
  }
}

