## Cool
load(file = "./inst/scenar_test_u=10.Rdata")
u10 <- compute_occurences(u10)

plot_diagram(u10, fill = "cnp", debug_mode = FALSE)
u10$run

#TODO: make a plot of multiple states
