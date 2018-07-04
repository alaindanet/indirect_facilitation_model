## Cool
load(file = "./inst/scenar_test_u=10.Rdata")
u10 <- compute_occurences(u10)

plot_diagram(u10, fill = "cnp", debug_mode = FALSE)
u10$run

#TODO: make a plot of multiple states

system.time(parallel::mclapply(1:24, function(i) Sys.sleep(1), mc.cores=24))
system.time(parallel::mclapply(1:24, function(i) Sys.sleep(1)))

run_scenarii_gradient()

u0 <- run_scenarii_gradient(
  gradient = list(g = seq(0, .3, length.out = 5), b = 1),
  model_spec = "two_facilitation_model",
  param = c(u = 0),
  time_seq = c(from = 0, to = 1, by = 1),
  solver_type = NULL
  )
avg_runs(u0)

cut_row = 2
u0
u0[["run"]] %>%
  dplyr::mutate(
    avg = parallel::mclapply(run, avg_runs, cut_row = cut_row)
    )
