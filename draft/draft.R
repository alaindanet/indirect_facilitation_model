library(simecol)
library(tidyverse)
library(magrittr)
devtools::load_all()
options(mc.cores = 1)

data(u0_5)


df_test <- filter(u0_5, g == .2, u == 5)
df_test2 <- df_test

df_test$run %<>% gather(species, rho, N, P) %>% 
  group_by(scenario, g, u, species) %>%
  nest() %>%
  mutate(
    group = map(data, identify_discontinuity, var = "rho")
    ) %>%
  unnest() %>%
  spread(species, rho)
high <- filter(df_test, scenario == "together")

unique(high$run$group)

ggplot(df_test$run,
  aes(x = b, y = rho, col = species,
    group = interaction(group, species))) +
  geom_line()

#TODO: modify plot bifurc to put interaction:
p <- plot_bifurcation(df_test2) +
  ggplot2::scale_colour_manual(
    values = c(N = "#BBCC33", P = "#99DDFF")
    ) + theme_alain()
p
debug(plot_bifurcation)


set.seed(123)
mod <- ca_two_facilitation_model()
solver(mod) <- myiteration3
times(mod) <- c(from = 0, to = 100, by = .5)
parms(mod)[c("u", "del", "gamma1", "b", "g", "f")] <-
  c(7.5, .4, 0.1, .9, .2, .9)
mod_run <- sim(mod)

load(file = "~/Documents/thesis/indirectfacilitation/inst/indirect_facilitation_and_dispersal_ca.RData")

gradient <- list(
  u = seq(7.5, 7.6, by = .1),
  del = seq(0.4, .41, by = .01)
  )

set.seed(123)
output <- run_scenarii_gradient(
  gradient = gradient,
  model_spec = "ca_two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1, b = .9, g =
.2, f = .9),
  time_seq = c(from = 0, to = 100, by = .5),
  set_tail = 300, nrep = 1
  )
save(output, file = "indirect_facilitation_and_dispersal_ca.RData")

#test <- output$run$run[[1]]
head(test)
head(output$run$run[[1]])


scenar_avg <- avg_runs(output, cut_row = 300)
rm(output)
save(scenar_avg, file = "indirect_facilitation_and_dispersal_ca_avg.RData")
rm(scenar_avg)
