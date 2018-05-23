# Effect of indirect and direct facilitation on coexistence and stability of a dryland ecosystem  

## TODO  

- [ ] Stability analysis:
    1. [x] Desert stability
	1. [x] Nurse seule
	2. [x] Protégé seule
	3. [x] Both species  
    2. [ ] Community stability
	1. [ ] Nurse seule
	2. [ ] Protégé seule
	3. [ ] Both species  
	4. [ ] P + 1 peu de N
	5. [ ] N + 1 peu de P
	6. [ ] 1 peu de P + 1 peu de N 
- [ ] Co-occurence analysis: f(b,g)
- [ ] Comparison pair-approx to mean field: f(b,g)
- [x] Test functions
    - [x] Test probability transitions
    - [ ] Test facilitate function
- [ ] Document functions
- [ ] Framework to build simulation diagrams
    - [ ] Check if simulations have reached stability
	- [x] Implementation of a custom solver ("lsodar")
	- [ ] Increase precision of the criterion [see here](https://stackoverflow.com/questions/5802592/dealing-with-very-small-numbers-in-r?rq=1) 
    - [x] Plot method 
	- [ ] print parameters of the simulations in the plot 
	    - [ ] Put an annotate method
	- [ ] Factorise plot_diagram (.gradient and .diagram): extract_state
	  function

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

