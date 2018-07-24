#!/usr/bin/env r


library(simecol)
library(tidyverse)
library(magrittr)
devtools::load_all()

options(mc.cores = 24)


set.seed(123)
gradient <- list(
  u = seq(0, 10, by = .1),
  del = seq(1, 0, by = -.01),
  g = c(0.1, .2)
  )
output <- run_scenarii_gradient(
  gradient = gradient,
  model_spec = "ca_two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1, b = .8),
  time_seq = c(from = 0, to = 6000, by = .5),
  set_tail = 10, nrep = 10
  )

save(output, file = "scenar_ca_cooccurence.Rdata")
scenar_avg <- avg_runs(output, cut_row = 10)
save(scenar_avg, file = "scenar_avg_ca_cooccurence.Rdata")
