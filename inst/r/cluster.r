#!/usr/bin/env r


library(simecol)
library(tidyverse)
library(magrittr)
devtools::load_all()

options(mc.cores = 24)

##########################
#  Complement to bifurc  #
##########################

gradient <- list(
  f = seq(0, 1, by = .01),
  del = seq(1, 0, by = -.01)
  )

output <- run_scenarii_gradient(
  gradient = gradient,
  model_spec = "two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1, g = .2, b = .9, u = 5),
  scenarii = init_scenarii(type = "together", ini_cover = .8),
  set_tail = 10
  )

options(mc.cores = 5)
scenar_avg <- avg_runs(output, cut_row = 10)
rm(output)
save(scenar_avg, file = "direct_facilitation_effect_pa_avg.RData")


#################
#  Cooccurence  #
#################

#Not run (around 120h with 20 cores)
options(mc.cores = 24)
set.seed(123)
output <- run_scenarii_gradient(
  gradient = gradient,
  model_spec = "ca_two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1, b = .9, g =
.2, u = 5),
  time_seq = c(from = 0, to = 10000, by = .5),
  set_tail = 300, nrep = 10
  )
save(output, file = "direct_facilitation_effect_ca_cooccurence.RData")

options(mc.cores = 5)
scenar_avg <- avg_runs(output, cut_row = 300)
rm(output)
save(scenar_avg, file = "direct_facilitation_effect_ca_avg.RData")
rm(scenar_avg)
