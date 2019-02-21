#!/usr/bin/env r


library(simecol)
library(tidyverse)
library(magrittr)
devtools::load_all()

options(mc.cores = 24)

###################################################
#  Stability along aridity and grazing gradients  #
###################################################

gradient <- list(
  g = seq(0, .5, by = .005),
  b = seq(1, 0.5, by = -.005),
  u = c(0, 5)
  )

output <- run_scenarii_gradient(
  gradient = gradient,
  model_spec = "two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1),
  scenarii = init_scenarii(type = c("bifurcation", "protegee_bifurcation"), ini_cover = .8, low_cover = 0.001),
  time_seq = c(from = 0, to = 10000, by = 1),
  set_tail = 10,
  solver_type = steady_state_4
  )

save(output, file = "bifurcation_diagram.RData")
options(mc.cores = 5)
scenar_avg <- avg_runs(output, cut_row = 10)
rm(output)
save(scenar_avg, file = "bifurcation_diagram_avg.RData")
