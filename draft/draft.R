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
