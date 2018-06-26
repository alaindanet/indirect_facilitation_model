# Effect of indirect and direct facilitation on coexistence and stability of a dryland ecosystem  

## TODO  

- [ ] Comparison pair-approx to mean field: f(b,g)
- [ ] Generalize run_scenarii_gradient and run_2d_gradient: f(b,g)
    - [ ] gradient argument: a list of **named** vectors
    - [ ] scenarii arg: NULL (default init values, scenario = "default"), a vector
      of valid starting values (scenario = "custom"), or a named list (provided by init_scenarii())
    - [ ] control for plot diagram, check how many param variables there is
- [ ] Update vignette
- [ ] Document functions
- [ ] Framework to build simulation diagrams
    - [ ] Check if simulations have reached stability
	- [x] Implementation of a custom solver ("lsodar")
	- [ ] Increase precision of the criterion [see here](https://stackoverflow.com/questions/5802592/dealing-with-very-small-numbers-in-r?rq=1) 
    - [ ] Put the number of gradient parameters that we want

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

