## Cool
load(file = "./inst/scenar_test_u=10.Rdata")


compute_states(u10, param = c("b", "g"), type = "double")

compute_states(u10, param = c(x = "b", y = "g"), type = "single")

states  <- c("coexistence", "nurse", "protégée", "extinct", "warning")
double_states <- expand.grid(states, states) %>%
      tidyr::unite(state) %>% unlist(.)
my_colours <- c(colors()[seq_along(double_states) - 1], "#C19A6B")

plot_diagram(u10,
  param = c(x = "b", y = "g"),
  possible_states = double_states,
  col_states = my_colours,
  debug_mode = FALSE, type = "double")
#TODO: make a plot of multiple states

init_scenarii(type = "together", model = two_facilitation_model())

