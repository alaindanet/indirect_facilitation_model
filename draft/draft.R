## Cool
library(simecol)
library(tidyverse)
library(magrittr)
devtools::load_all()
options(mc.cores = 1)

load(file = "inst/scenar_avg_bifurc_u=0_5_gamma1_.1.Rdata")

states <- compute_states(scenar_avg, type = "double")
plot_diagram(states) + theme_get()

df <- filter(scenar_avg, g == .2, u == 5)

p <- plot_bifurcation(df) +
  ggplot2::scale_colour_manual(
    values = c(N = "#BBCC33", P = "#99DDFF")
    ) + theme_alain()
p

df_high <- filter(df, scenario == "together")

ind <- abs(
  as.numeric(df_high$run$N[-1]) - as.numeric(df_high$run$N[-nrow(df_high$run)])
  ) >= 0.1
splitAt <- function(x, pos) split(x, cumsum(seq_along(x) %in% (pos+1)))
l1 <- splitAt(as.numeric(df_high$run$N), which(ind))
names(l1) <- 1:length(l1)
l2 <- lapply(seq_along(l1), 
             function(y, n, i) {
                                 as.numeric(rep(n[[i]], length(y[[i]]))) 
                               }, y=l1, n=names(l1))
unlist(l2)

myxylab <- function (...) {
  dots <- pryr::named_dots(...)
  dots <- unlist(dots)

  lab_list <- list(
    mpg = paste("Miles/(US) gallon"),
    qsec = expression(paste("1/4 mile time (", tau, ")")),
    disp = paste("Displacement (cu.in.)"),
    wt = paste("Weight (1000 lbs) ")
    )
  
  lab_used <- lab_list[dots]
  names(lab_used) <- names(dots)

  labs(
    x = lab_used["x"][[1]],
    y = lab_used["y"][[1]]
    )

}

qplot(mpg, wt, data = mtcars) +
  labs(
    x = paste("Miles/(US) gallon"),
    y = paste("Weight (1000 lbs)")
    )

qplot(mpg, wt, data = mtcars) +
  myxylab(x = "mpg", y = "wt")

