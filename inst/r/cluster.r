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
  u = seq(0, 10, by = .1),
  del = seq(1, 0, by = -.01)
  )

output <- run_scenarii_gradient(
  gradient = gradient,
  model_spec = "two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1, g = .2, b = .9),
  scenarii = init_scenarii(type = "together", ini_cover = .8),
  set_tail = 10
  )

options(mc.cores = 5)
scenar_avg <- avg_runs(output, cut_row = 10)
rm(output)
save(scenar_avg, file = "clustering_pa_avg.RData")


#################
#  Cooccurence  #
#################

#Not run (around 120h with 20 cores)
set.seed(123)
output <- run_scenarii_gradient(
  gradient = gradient,
  model_spec = "ca_two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1, b = .9, g = .2),
  time_seq = c(from = 0, to = 10000, by = .5),
  set_tail = 300, nrep = 10
  )
save(output, file = "scenar_ca_cooccurence.Rdata")

options(mc.cores = 5)
scenar_avg <- avg_runs(output, cut_row = 300)
rm(output)
save(scenar_avg, file = "clustering_ca_avg.Rdata")
rm(scenar_avg)
