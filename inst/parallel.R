#!/usr/bin/env r

################################################################################
#                       R script for parrallel computing                       #
################################################################################

# Library loading
needed_library <- c("tidyverse", "multidplyr", "simecol", "magrittr", "devtools")
sapply(needed_library, function(x){ library(x, character.only = TRUE) })
# load package:
devtools::load_all()

################################################################################
# Set parameters
g_gradient <- seq(0, 0.3, length.out = 100)
b_gradient <- seq(1, 0.5, length.out = 100)
cores <- 20 

################################################################################
# Simulation

# u0
u0 <- run_scenarii_gradient( y = "g",
  gradienty = g_gradient,
  x = "b",
  gradientx = b_gradient,
  model_spec = "two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1, u = 0),
  time_seq = c(from = 0, to = 3000, by = 1),
  nb_cores = cores,
  solver_type = NULL,
  scenarii = init_scenarii(type = "bifurcation", ini_cover = .8, low_cover = .01)
  ) %>% avg_runs(.)
save(u0, file = "scenar_bifurc_u=0.Rdata")

# u5
u5 <- run_scenarii_gradient( y = "g",
  gradienty = g_gradient,
  x = "b",
  gradientx = b_gradient,
  model_spec = "two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1, u = 5),
  time_seq = c(from = 0, to = 3000, by = 1),
  nb_cores = cores,
  solver_type = NULL,
  scenarii = init_scenarii(type = "bifurcation", ini_cover = .8, low_cover = .01)
  ) %>% avg_runs(.)
save(u5, file = "scenar_bifurc_u=5.Rdata")

