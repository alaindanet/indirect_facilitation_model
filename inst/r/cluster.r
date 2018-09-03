#!/usr/bin/env r


library(simecol)
library(tidyverse)
library(magrittr)
devtools::load_all()

options(mc.cores = 20)

##########################
#  Complement to bifurc  #
##########################

gradient <- list(
  u = 10,
  b = seq(1, 0.5, by = -0.005),
  g = seq(0, .5, 0.005)
  )

output <- run_scenarii_gradient(
  gradient = gradient,
  model_spec = "two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1),
  scenarii = init_scenarii(type = "bifurcation", ini_cover = .8, low_cover = .1),
  set_tail = 10
  )

options(mc.cores = 5)
scenar_avg <- avg_runs(output, cut_row = 10)
rm(output)
save(scenar_avg, file = "scenar_birfuc_u_10.RData")
