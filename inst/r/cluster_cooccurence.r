#!/usr/bin/env r


library(simecol)
library(tidyverse)
library(magrittr)
devtools::load_all()


#########################
#  Direct facilitation  #
#########################

##Not run (around 120h with 20 cores)
#options(mc.cores = 24)
## Define parameter gradient
#gradient <- list(
#  f = seq(0, 1, by = .01),
#  del = seq(0, 1, by = .01)
#  )
#
#set.seed(123)
#output <- run_scenarii_gradient(
#  gradient = gradient,
#  model_spec = "ca_two_facilitation_model",
#  param = c(protection_type = list("first_protect"), gamma1 = 0.1, b = .9, g =
#.2, u = 5),
#  time_seq = c(from = 0, to = 10000, by = .5),
#  set_tail = 300, nrep = 10
#  )
#save(output, file = "direct_facilitation_and_dispersal_ca.RData")
#
#options(mc.cores = 5)
#scenar_avg <- avg_runs(output, cut_row = 300)
#rm(output)
#save(scenar_avg, file = "direct_facilitation_and_dispersal_ca_avg.RData")
#rm(list = ls())


###########################
#  Indirect facilitation  #
###########################

#Not run (around 120h with 20 cores)
options(mc.cores = 24)
# Define parameter gradient
gradient <- list(
  u = seq(0, 10, by = .1),
  del = seq(0, 1, by = .01)
  )

set.seed(123)
output <- run_scenarii_gradient(
  gradient = gradient,
  model_spec = "ca_two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1, b = .9, g =
.2, f = .9),
  time_seq = c(from = 0, to = 10000, by = .5),
  set_tail = 300, nrep = 10
  )
save(output, file = "indirect_facilitation_and_dispersal_ca.RData")

options(mc.cores = 5)
scenar_avg <- avg_runs(output, cut_row = 300)
rm(output)
save(scenar_avg, file = "indirect_facilitation_and_dispersal_ca_avg.RData")
rm(scenar_avg)
