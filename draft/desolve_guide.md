
# When desolve send warnings

## Increase timestep 

You can do it by decreasing the `by` argument of the `time` argument:

```r
c(from = 0, to = 100, by = 0.1)
```

## Debug with the "rk4" solver

(Reference)[http://r.789695.n4.nabble.com/Problems-with-the-deSolve-package-td4717598.html].

The rk4 solver allow to run simulations with a constant bandwidth, so you can
see what is going wrong. 

```r
solver(model) <- "rk4" 
```

