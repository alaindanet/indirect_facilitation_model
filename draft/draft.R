## Cool
load(file = "./inst/scenar_test_u=10.Rdata")
u10 <- compute_occurences(u10)

plot_diagram(u10, fill = "cnp", debug_mode = FALSE)
u10$run

#TODO: make a plot of multiple states


run_scenarii_gradient()

run_scenarii_gradient(
  gradient = list(g = 1, b = 1),
  model_spec = "two_facilitation_model",
  time_seq = c(from = 0, to = 1, by = 1),
  solver_type = NULL
  )
split(t(param_combination)[[1]], f = names(param_combination))

switch_list <- function (l) {
  apply(l, 1, function(x){})
  
}


df2list(param_combination)


Map(dummy, comb[["inits"]],
  param_combination,
  MoreArgs = list(z = "two_facilitation_model"))
param_test <- list(c(b = 1, c = 0.2))

Map(dummy, comb[["inits"]][1],
  param_test,
  MoreArgs = list(z = "two_facilitation_model"))

