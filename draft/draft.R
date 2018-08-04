## Cool
library(simecol)
library(tidyverse)
library(magrittr)
options(mc.cores = 1)


scenar_occurences <- compute_occurences(scenar_avg)$run %>%
  gather(var, clustering, cnp, cnn, cpp, cveg, cveg2) %>%
  group_by(g, b, var) %>%
  summarise(clustering = mean(clustering, na.rm = TRUE))
scenar_occurences

load(file = "inst/scenar_avg_ca_cooccurence.Rdata")
scenar_avg$run <- tibble::as.tibble(scenar_avg$run)

occurences <- scenar_avg
occurences$run <- filter(scenar_avg, P > 0.01 | P == 0, N > 0.01 | N == 0)  %>%
  compute_occurences(.) %>%
  .$run %>%
  gather(clustering, value, cnp:cveg) %>%
  group_by(u, del, g, clustering) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  spread(clustering, value)

# Add log10 scale to the plot
# Add color scale with diverging palette 
plot_fig3(occurences)

## Check big cooccurences
big_co_occcurences <- filter(occurences, cpp > 100)

gradient_test <- sapply(big_co_occcurences$gradient, function(x) x[1])
gradient_test["del"] <- big_co_occcurences$gradient$del[4]
gradient_test["del"] <- .1

output <- run_scenarii_gradient(
  gradient = gradient_test,
  model_spec = big_co_occcurences$model,
  param = big_co_occcurences$param,
  time_seq = c(from = 0, to = 10000, by = .5),
  set_tail = 300, nrep = NULL
  )
avg_runs(output, cut_row = 300)


mod <- ca_two_facilitation_model()
solver(mod) <- "ca_solver"
times(mod)["to"] <- c(to = 20000)
parms(mod)["g"] <- .1
set.seed(123)
mod_run <- sim(mod)
#plot(mod_run)

g1 <- out(mod_run) %>%
  gather(state, rho, N, P, E, D) %>%
  ggplot(., aes(x = time, y = rho, color = state)) +
  geom_line() + ylim(0, 1)
g1

test <- out(mod_run)
test[204,]$P == 0
