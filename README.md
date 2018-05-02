# Effect of indirect and direct facilitation on coexistence and stability of a dryland ecosystem  

## TODO  

- [x] Modularize:  
    - [x] make easy the specification of the associative protection
- [x] Test functions
    - [ ] Compare new simulations to master's one
    - [x] Test probability transitions
    - [ ] Test facilitate function
    - [ ] Test that the parameters are well spelled       
- [ ] Document functions
- [ ] Framework to build simulation diagrams
    - [ ] Add a class "gradient" with appropriate methods 
	- [x] Build `run_bifurcation` (`run_2d_gradient` wrapper)
	- [x] Build `plot_np` method for bifurcation class
    - [ ] Check if simulations have reached stability
	- [x] Implementation of a custom solver ("lsodar")
	- [ ] Increase precision of the criterion [see here](https://stackoverflow.com/questions/5802592/dealing-with-very-small-numbers-in-r?rq=1) 
    - [x] Define a plot method 
	- [ ] print parameters of the simulations in the plot 
	    - [ ] Put an annotate method
	    - [x] Keep the parameters values in the simulation Framework

## Objectives

Develop a theoretical model of two species in competition and add an interaction
of indirect facilitation. A nurse species is protected against grazing can protect the saplings of a protégée species by associative protection. We make the simple assumption of a growth-defence trade-off.

The cells can take 3 states: occupied by a nurse, a protégée or be empty.

## Informations about the model

- 3 state model:
- Grazing protection specification

## Workflow (development)

### Clone the repo

```
git clone https://github.com/alaindanet/indirect_facilitation_model.git my_dir_name

cd my_dir_name
```

### Go

- Load the packages:

```
devtools::load_all()

library(simecol)
library(ggplot2)
library(tibble)
library(magrittr)
```

- Load the model and specify the parameters:

```
mod <- indirect_facilitation_model()
mod

# Show parameters:
parms(mod)

# Tweak parameters:
parms(mod)["g"] <- 0
parms(mod)["gamma1"] <- 0.20
times(mod) <- c(from = 0, to = 1000, by = 1)
```

- Launch a simulation and plot the result:

```
mod_run <- sim(mod)
plotnp(mod_run)
```

