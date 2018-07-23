## Cool
library('simecol')


mod <- ca_two_facilitation_model()
solver(mod) <- "myiteration"
times(mod)["to"] <- c(to = 1000)
parms(mod)["g"]  <- c(.1)
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
parms(mod)["g"]  <- c(.1)


u0 <- run_scenarii_gradient(
    gradient = list(g = c(0, 0.1), b = c(.5, .8)),
    model_spec = "ca_two_facilitation_model",
    time_seq = c(from = 0, to = 1, by = 1),
    solver_type = NULL, nrep = 3
  )
avg_runs(u0)
