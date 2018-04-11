library(simecol)
library(ggplot2)
library(tibble)
library(magrittr)

devtools::document()

upca <- upca_model()
parms(upca)
equations(upca)$f <- equations(upca)$f2
test <- sim(upca)
plotupca(test)

##   

# Test
mod <- indirect_facilitation_model()
mod
times(mod) <- c(from = 0, to = 300, by = .3)
equations(mod)$asymp_cost <- equations(mod)$asymp_cost
mod_run <- sim(mod)

parms(mod)["g"] <- 0
parms(mod)["gamma1"] <- 0.20
times(mod) <- c(from = 0, to = 1000, by = 1)
mod_run <- sim(mod)
tail(out(mod_run), 40)
plotnp(mod_run)

# Run_gradient
g_gradient <- seq(0, 0.2, length.out = 10)
gamma_gradient <- seq(0, 0.2, length.out = 10)

spec_param <- parms(mod)["b"] <- .8

gradient_2d <- run_2d_gradient(gradienty = g_gradient,
  time_seq = c(from = 0, to = 1000, by = 1), nb_cores = 4
  )

gradient_2d

averaged_runs <- gradient_2d %>%
  dplyr::mutate(avg = purrr::map(runs, avg_runs, cut_row = 10)) %>%
  tidyr::unnest(avg)

plot_diagram(averaged_runs)
ggsave("./inst/figs/diag_cost_grazing.pdf")
plot_diagram(averaged_runs, debug_mode = TRUE)
# Not work bc many values by gradient value 
plotnp_gradient(averaged_runs)

# filter
test <- dplyr::filter(averaged_runs, g == 0 & gamma1 == 0.05)
with(test, {
  def_state(N,P,status)
  })

test2 <- dplyr::filter(gradient_2d, g == 0 & gamma1 == 0.2)
tail(test2[1,]$runs[[1]])
