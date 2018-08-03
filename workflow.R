library(simecol)
library(tidyverse)
library(magrittr)
options(mc.cores = 3)

devtools::document()
#devtools::use_vignette("three_states_model")

mod <- indirect_facilitation_model()
init(mod)
parms(mod)["g"] <- 0.13
parms(mod)["gamma1"] <- 0.1
parms(mod)["protection_type"] <- list("first_protect")
parms(mod)["protection_type"] <- list("linear")
times(mod) <- c(from = 0, to = 3000, by = 1)
mod_run <- sim(mod)
tail(out(mod_run), 10)
plotnp(mod_run)

################################
#  Grazing and gamma gradient  #
################################

# Run_gradient g and gamma1
gradient <- list(
  g = seq(0, 0.3, length.out = 10),
  gamma1 = seq(0, 0.3, length.out = 10)
  )
gradient_2d <- run_scenarii_gradient(
  gradient = gradient,
  param = c(b = 0.8, u = 0, tau_n = 0, protection_type = list("linear")),
  time_seq = c(from = 0, to = 3000, by = 1),
  solver_type = NULL,
  scenarii = init_scenarii(type = "together")
  )

gradient_2d
averaged_runs <- avg_runs(gradient_2d, cut_row = 10)

plot_diagram(averaged_runs, param = c(x = "gamma1", y = "g"), debug_mode = FALSE)
#ggsave("diag_g_gamma_linear_tau_n=0.png")
plot_diagram(averaged_runs, debug_mode = TRUE)

##################################
#  Grazing and aridity gradient  #
##################################

gradient <- list(
  g = seq(0, 0.3, length.out = 10),
  b = seq(1, 0.5, length.out = 10)
  )
gradient_2d <- run_scenarii_gradient(
  gradient = gradient,
  param = c(protection_type = list("first_protect"), gamma1 = 0.1, u = 0),
  time_seq = c(from = 0, to = 3000, by = 1),
  scenarii = init_scenarii(type = "together")
  )

averaged_runs <- avg_runs(gradient_2d, cut_row = 10)
plot_diagram(averaged_runs, param = c(x = "b", y = "g"))
ggsave("./inst/figs/diag_aridity_grazing_first_protect_u=0.png")

###############################
#  Bifurcation state diagram  #
###############################

bifurc <- run_bifurcation(
  gradientx = seq(0.2, 1, by = 0.1),
  gradienty = c(0.4, 0.01),
  param = c(
    g = .15, gamma1 = .1,
    protection_type = list("first_protect"), #list("first_protect")
    u = 0 
    )
  )

averaged_runs <- avg_runs(bifurc, cut_row = 1)
plotnp(averaged_runs, alpha = 0.65)

################################################################################
#                              Four states model                               #
################################################################################

mod <- two_facilitation_model()
init(mod)
names(init(mod))
parms(mod)["g"] <- 0.3
parms(mod)["gamma1"] <- 0.08
parms(mod)["protection_type"] <- list("first_protect")
parms(mod)["protection_type"] <- list("linear")
times(mod) <- c(from = 0, to = 1000, by = 1)
parms(mod)["u"] <- 30
mod_run <- sim(mod)
plotnp(mod_run)

#############################################
#  Grazing and grazing protection strength  #
#############################################

gradient <- list(
  g = seq(0, 0.3, length.out = 10),
  u = seq(0, 20, length.out = 10)
  )

gradient_2d <- run_scenarii_gradient(
  gradient = gradient,
  model_spec = "two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1),
  time_seq = c(from = 0, to = 3000, by = 1),
  scenarii = init_scenarii(type = "together")
  )

averaged_runs <- avg_runs(gradient_2d, cut_row = 10)
plot_diagram(averaged_runs, param = c(x = "u", y = "g"))
ggsave("inst/figs/four_states/diag_u_grazing_first_protect_gamma1=.1.png")

##################################
#  Grazing and aridity gradient  #
##################################

gradient <- list(
  g = seq(0, 0.3, length.out = 30),
  b = seq(1, .5, length.out = 30)
  )
gradient_2d <- run_scenarii_gradient(
  gradient = gradient,
  model_spec = "two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1, u = 5),
  time_seq = c(from = 0, to = 3000, by = 1),
  scenarii = init_scenarii(type = "together")
  )

averaged_runs <- avg_runs(gradient_2d, cut_row = 10)
plot_diagram(averaged_runs, param = c(x = "b", y = "g"))
#save(gradient_2d, file = "diag_aridity_grazing_first_protect_u=5.RData" )
ggsave("inst/figs/four_states/diag_aridity_grazing_first_protect_u=5.pdf",
  scale = .8)

###############################
#  Bifurcation state diagram  #
###############################

bifurc <- run_bifurcation(
  gradientx = seq(0.2, 1, length.out = 30),
  gradienty = c(0.4, 0.01),
  model_spec = "two_facilitation_model",
  time_seq = c(from = 0, to = 3000, by = 1),
  param = c(
    g = .3, gamma1 = .1,
    protection_type = list("first_protect"), #list("first_protect")
    u = 5
    )
  )

averaged_runs <- avg_runs(bifurc, cut_row = 1)
plotnp(averaged_runs, alpha = 0.65)
ggsave("inst/figs/four_states/bifurc_first_protect_u=0_g=.25_gamma1=.1.png")

##############
#  Scenarii  #
##############

gradient <- list(
  g = seq(0, 0.3, length.out = 30),
  b = seq(1, .5, length.out = 30)
  u = c(0, 5)
  )
output <- run_scenarii_gradient(
  gradient = gradient,
  model_spec = "two_facilitation_model",
  param = c(protection_type = list("first_protect"), gamma1 = 0.1),
  time_seq = c(from = 0, to = 3000, by = 1),
  scenarii = init_scenarii(type = "bifurcation", ini_cover = .8, low_cover = .01)
  )

save(output, file = "./inst/scenar_bifurc_u=0_5_gamma1_.1.Rdata")
scenar_avg <- avg_runs(output)
save(scenar_avg, file = "./inst/scenar_avg_bifurc_u=0_5_gamma1_.1.Rdata")

g <- plot_diagram(scenar_avg, param = c(x = "b", y = "g", type = "scenario"),
  debug_mode = FALSE)
g

test <- compute_states(scenar_avg, param = c(x = "b", y = "g", type = "scenario"))
test
test$run %>%
  spread(scenario, state)

plotnp(scenar_avg)
ggsave("inst/figs/four_states/scenar_bifurc_first_protect_u=0_test.png")

names(output)

test <- avg_runs(output$u10)
plot_diagram(test)

##############
#  Figure 1  #
##############
library("stringr")

load(file = "./inst/scenar_bifurc_u=0_10_15_gamma1_.1.Rdata")
test <- lapply(output, FUN = avg_runs)
rm(output)
test2 <- lapply(test, FUN = compute_states, param = c(x = "b", y = "g", type = "scenario"))

test3 <- list()
for (i in seq_along(names(test2))) {
  test3[[i]] <- test2[[names(test2)[i]]] %>%
    dplyr::mutate(
      u = as.numeric(substr(names(test2)[i], 2, 10)),
      state = stringr::str_replace(state, "warning", "extinct")
      )
}
allu <- bind_rows(test3)

# À gauche 
u0 <- allu %>%
  filter(scenario == "together", u == 0, b >= 0.5)
class(u0) <- c("tibble", "data.frame", "states")
plot_diagram(u0, param = c(x = "b", y = "g"))
ggsave("inst/figs/four_states/diag_u=0.pdf", width = 7, height = 5, units =
  "cm")

# À droite 
load("diag_aridity_grazing_first_protect_u=5.RData") # gradient_2d

u5 <- avg_runs(gradient_2d)

compute_states(u5, param = c(x = "b", y = "g"))$run

states_u5 <- compute_states(u5, param = c(x = "b", y = "g")) %>%
  filter(b >= 0.5)
class(states_u5) <- c("tibble", "data.frame", "states")
plot_diagram(states_u5, param = c(x = "b", y = "g"))
ggsave("inst/figs/four_states/diag_u=5.pdf", width = 7, height = 5, units = "cm")


# Bonus
u10 <- allu %>%
  filter(scenario == "together", u == 10)
class(u10) <- c("tibble", "data.frame", "states")
plot_diagram(u10, param = c(x = "b", y = "g"))
ggsave("inst/figs/four_states/diag_u=10.pdf", scale = .8)

##############
#  Figure 2  #
##############

#u_grad <- c(0, 5)
u_grad <- c(0)
#g_grad <- c(0, .05, .1, .2, .3)
g_grad <- c(.09)

for(u in seq_along(u_grad)){
  for(g in seq_along(g_grad)) {

    bifurc <- run_bifurcation(
      gradientx = seq(0.2, 1, length.out = 30),
      gradienty = c(0.4, 0.01),
      model_spec = "two_facilitation_model",
      time_seq = c(from = 0, to = 3000, by = 1),
      param = c(
	g = g_grad[g], gamma1 = .1,
	protection_type = list("first_protect"),
	u = u_grad[u]
	)
      )

    averaged_runs <- avg_runs(bifurc, cut_row = 1)
    plotnp(averaged_runs, alpha = 0.65)
    plot_name <- paste("inst/figs/four_states/bifurc_first_protect_u=",
      u_grad[u], "_g=", g_grad[g], "_gamma1=.1.pdf", sep = "")
    ggsave(plot_name, width = 7, height = 5, units = "cm")
  }
}

#######################
#  Figure 2 improved  #
#######################

load(file = "./inst/scenar_test_u=10.Rdata")
double_states <- c("coexistence", "nurse", "protegee", "desert",
  "protegee_desert", "nurse_desert", "coexistence_desert", "protegee_nurse",
  "unkown")
my_colours <- c(coexistence = "orange", nurse = "green", protegee = "black",
  desert = "#C19A6B", protegee_desert = "gray60",
  nurse_desert = "green2", coexistence_desert = "orange2", protegee_nurse =
    "darkgreen", unkown = "white")

plot_diagram(u0,
  param = c(x = "b", y = "g"),
  possible_states = double_states,
  col_states = my_colours,
  debug_mode = FALSE, type = "double")
ggsave(filename = "./inst/figs/four_states/diag_bistab_u=10_gamma1=.1.pdf")


######################
#  Figure 2 new sim  #
######################
load(file = "./inst/scenar_avg_bifurc_u=0_5_gamma1_.1.Rdata")

states <- compute_states(scenar_avg, type = "double")

plot_fig2(states)
ggsave("./inst/figs/figure2.pdf", device = cairo_pdf)


###########################################
#  Plot species densities and clustering  #
###########################################

plot_diagram(scenar_avg, debug_mode = FALSE, fill = "P")

occurenres <- compute_occurences(scenar_avg)
summary(occurenres$run$c_veg3)
plot_diagram(occurenres, fill = "cnp")
plot_diagram(occurenres, fill = "cpp")
plot_diagram(occurenres, fill = "cnn")
plot_diagram(occurenres, fill = "c_veg3")

######################
#  Bifucation plots  #
######################


extraction <- filter(scenar_avg, g %in% c(0, 0.1, 0.20, .3))

p <- plotnp(extraction, b, threshold = 10^-3, debug_mode = FALSE, N, P)
p + facet_grid(vars(g), vars(u)) 

################################################################################
#                              Cellular automata                               #
################################################################################


mod <- ca_two_facilitation_model()
solver(mod) <- "myiteration"
times(mod)["to"] <- c(to = 6000)
parms(mod)
mod_run <- sim(mod)
#plot(mod_run)

g1 <- out(mod_run) %>%
  gather(state, rho, -times) %>%
  ggplot(., aes(x = times, y = rho, color = state)) +
  geom_line() + ylim(0, 1)
g1

#################
#  Cooccurence  #
#################

load(file = "inst/scenar_avg_ca_cooccurence.Rdata")
scenar_avg$run <- tibble::as.tibble(scenar_avg$run)

# Average by replicate
occurences <- scenar_avg
occurences$run <- compute_occurences(scenar_avg) %>%
  .$run %>%
  gather(clustering, value, cnp:cveg) %>%
  group_by(u, del, g, clustering) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  spread(clustering, value)

plot_list <- list()
clustering <- c("cnp", "cpp", "cnn", "cveg")
for (i in seq_along(clustering)) {
  plot_list[[i]] <- plot_diagram(occurences, param = c(x = "del", y = "u"), fill = clustering[i])
}

cowplot::plot_grid(plotlist = plot_list, labels = "auto")
