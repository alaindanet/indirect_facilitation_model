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

sapply(as.list(out(mod_run)[, c("nurse", "protegee")]), moments::skewness)

while (x & i < 10 & i < 30 & y) {
  i <- i + 1
  print(i)
  if(i == 5){ x <- FALSE; y  <- FALSE}
}
 

mod <- ca_two_facilitation_model()
times(mod)["to"] <- c(to = 6000)
parms(mod)["g"]  <- c(.1)
library(microbenchmark)


