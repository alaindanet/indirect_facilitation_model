## Cool
load(file = "./inst/scenar_test_u=10.Rdata")
u10 <- compute_occurences(u10)

plot_diagram(u10, fill = "cnp", debug_mode = FALSE)
u10$run

#TODO: make a plot of multiple states

system.time(parallel::mclapply(1:3, function(i) Sys.sleep(1), mc.cores = 3))
system.time(parallel::mclapply(1:3, function(i) Sys.sleep(1)))

run_scenarii_gradient()

u0 <- run_scenarii_gradient(
  gradient = list(
    g = 0,
    b = seq(1, .5, length.out = 10)
    ),
  model_spec = "two_facilitation_model",
  param = c(u = 0),
  time_seq = c(from = 0, to = 3000, by = 1),
  scenarii = init_scenarii(type = "bifurcation")
  )
u0_avg <- avg_runs(u0)

plotnp(u0_avg, b, threshold = 10^-3, debug_mode = FALSE, N, P)
filter(u0_avg$run, N > 0)
filter(u0_avg, N > 0)

names(u0_avg$inits[[1]])
