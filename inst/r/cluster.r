#!/usr/bin/env r


library(simecol)
library(tidyverse)
library(magrittr)
devtools::load_all()

options(mc.cores = 24)


gradient <- list(
  g = seq(0, 0.5, length.out = 100) %>% round(., 3),
  b = seq(1, .5, length.out = 100) %>% round(., 3),
  u = c(0, 5)
  )
output <- run_scenarii_gradient(
  gradient = gradient,
  model_spec = "two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1),
  time_seq = c(from = 0, to = 5000, by = 1),
  scenarii = init_scenarii(type = "bifurcation", ini_cover = .8, low_cover = .01),
  set_tail = 10
  )

save(output, file = "./inst/scenar_bifurc_u=0_5_gamma1_.1.Rdata")
scenar_avg <- avg_runs(output, cut_row = 10)
save(scenar_avg, file = "./inst/scenar_avg_bifurc_u=0_5_gamma1_.1.Rdata")
