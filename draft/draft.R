## Cool
library(simecol)
library(tidyverse)
library(magrittr)
options(mc.cores = 3)


mod <- ca_two_facilitation_model()
times(mod)["to"] <- c(to = 1000)
parms(mod)[c("g", "b", "del")]  <- c(.1, .4, 1)
mod_run <- sim(mod)
#plot(mod_run)

test <- out(mod_run)
sapply(test, moments::skewness)

g1 <- as.tibble(out(mod_run)) %>%
  mutate(
    cut_time = c(0,sapply(seq(n() / 300), function(x) {
      rep(x, 300)})),
    cut_time = as.integer(cut_time)
    ) %>%
  gather(state, rho, nurse, protegee) %>%
  filter(cut_time != 0) %>%
  group_by(cut_time, state) %>%
  summarise(skewness = moments::skewness(rho))
ggplot(g1, aes(x = cut_time, y = skewness, fill = state)) +
  geom_point()
g1

mod <- ca_two_facilitation_model()
times(mod)["to"] <- c(to = 6000)
parms(mod)[c("g", "b")]  <- c(.1, .5)
sim(mod)

scenar_occurences <- compute_occurences(scenar_avg)$run %>%
  gather(var, clustering, cnp, cnn, cpp, cveg, cveg2) %>%
  group_by(g, b, var) %>%
  summarise(clustering = mean(clustering, na.rm = TRUE))
scenar_occurences

load(file = "inst/scenar_avg_ca_cooccurence.Rdata")
scenar_avg$run <- tibble::as.tibble(scenar_avg$run)

occurences <- scenar_avg
occurences$run <- compute_occurences(scenar_avg) %>%
  .$run %>%
  gather(clustering, value, cnp:cveg) %>%
  group_by(u, del, g, clustering) %>%
  summarise(value = mean(value, na.rm = TRUE) %>% log10(.)) %>%
  spread(clustering, value)

plot_fig3(occurences)

clustering <- occurences
clustering$run %<>% gather(var, c, cnn, cnp, cpp, cveg)
plot_diagram(clustering, param = c(x = "del", y = "u"), fill = "c") +
  scale_fill_gradient2(low = scales::muted("blue"), mid = "white",
    high = scales::muted("red")) +
  facet_grid(vars(g), vars(var), labeller = labeller(g = as_labeller(g_appender), var = cxx))

filter(occurences, cpp > 3)

#########################
#  Averaging densities  #
#########################
occurences <- compute_occurences(scenar_avg)
occurences$run %<>%
  gather(rho, value, N:qveg) %>%
  group_by(u, del, g, rho) %>%
  summarise(value = mean(value, na.rm = TRUE)%>% log10(.)) %>% 
  spread(rho, value)
occurences$run <- compute_occurences(occurences) %>%
  .$run %>%
  gather(clustering, value, N:cveg) %>%
  group_by(u, del, g, clustering) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  spread(clustering, value)

clustering <- occurences
clustering$run %<>% gather(var, c, cnn, cnp, cpp, cveg)
plot_diagram(clustering, param = c(x = "del", y = "u"), fill = "c") +
  scale_fill_gradient2(low = scales::muted("blue"), mid = "white", high = scales::muted("red")) +
  facet_grid(vars(g), vars(var))
