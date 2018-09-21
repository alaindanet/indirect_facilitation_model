## Cool
library(simecol)
library(tidyverse)
library(magrittr)
options(mc.cores = 1)

load(file = "inst/scenar_avg_bifurc_u=0_5_gamma1_.1.Rdata")

load(file = "./inst/scenar_avg_bifurc_u=0_5_gamma1_.1.Rdata")

states <- compute_states(scenar_avg, type = "double")
plot_diagram(states) + theme_get()

