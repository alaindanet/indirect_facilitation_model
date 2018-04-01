
u <- .4
k <- 2
curve(k * x ^ u, 0, .2)

## Asymptote
tau <- 20
asymp <- 1
curve(1 - exp(-x * tau), 0, .2, 
  xlab = "gamma", 
  ylab = "n")
