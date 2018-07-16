## Cool
load(file = "./inst/scenar_bifurc_u=0_5_gamma1_.1.Rdata")

avg_runs(output, cut_row = 5)

u0 <- filter(scenar_avg, u == 0)
u5 <- filter(scenar_avg, u == 5)
names(u0$gradient)

## TODO: check the computation of states
states <- compute_states(scenar_avg, type = "double")
unique(states$run$state)
which(is.na(states$run$state), arr.ind = TRUE)

plot_diagram(states, type = "double_states", debug_mode = FALSE) +
    ggplot2::scale_fill_manual(
      values = color_states()
      ) + facet_grid(cols = vars(u))


