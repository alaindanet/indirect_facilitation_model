# Effect of indirect and direct facilitation on coexistence and stability of a dryland ecosystem  

## TODO  

- [ ] Modularize:  
    - [ ] make easy to the specification of the associative protection
- [x] Test functions
    - [x] Test probability transitions
    - [ ] Test  
- [ ] Document functions
- [ ] Framework to build simulation diagrams
    - [x] Run simulations over a range of parameters  
    - [x] Extract the average of the last timesteps
    - [x] Check if simulations become weird (return a WARNING ?)
    - [x] Ajouter un seuil d'extinction
    - [ ] Check if simulations have reached stability
	- [x] Implementation of a custom solver ("lsodar")
	- [ ] Increase precision of the criterion [see here](https://stackoverflow.com/questions/5802592/dealing-with-very-small-numbers-in-r?rq=1) 
    - [x] Define a plot method 
	- [ ] print parameters of the simulations in the plot 
	- [x] Generalize the 2d gradient 
	- [x] Generalize the plot method 

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
