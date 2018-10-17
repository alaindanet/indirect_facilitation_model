# Effect of indirect and direct facilitation on coexistence and stability of a dryland ecosystem  

## TODO  

### Paper 

- Fig 1: methods
- Fig 2: multi-states + bifurcation details
- Fig 3: cellular automata, spatial coexistence Cnp $f(\gamma, u)$ 

### Implementation

- [ ] sim_multi():
    - [ ] Save a simecol object
    - [ ] Save a time argument
    - [ ] Save init values
    - [ ] Save baseline parameters 
    - [ ] Save param_combination
    - [ ] param_combination: matrix of parameter combination
    - [ ] be able to sim from a sim_multi object
- [ ] Generalize run_scenarii_gradient and run_2d_gradient: f(b,g)
    - [ ] check_param(model_spec, names(gradient)))
    - [ ] Add check for argument conformity (check_param, check_inits)
      of valid starting values (scenario = "custom"), or a named list (provided by init_scenarii())
    - [ ] Add methods: init.scenarii(), param.scenarii() (S4 methods)
      https://stackoverflow.com/questions/12100856/combining-s4-and-s3-methods-in-a-single-function
- [ ] plot: separate lines (create groups: e.g: together & < threshold, together & >= threshold, low_together ...)
    - [ ] See [here](https://stackoverflow.com/a/23863893/5968131) 
- [x] Compare the effect of dispersal, facilitation strength, paturâge, aridité on clustering  
- [ ] Document functions
- [x] Check if simulations have reached stability
    - [x] Implementation of a custom solver ("lsodar")
    - [x] Test the implementation

## Objectives

Develop a theoretical model of two species in competition and add an interaction
of indirect facilitation. A nurse species is protected against grazing can protect the saplings of a protégée species by associative protection. We make the simple assumption of a growth-defence trade-off.


## Informations about the model

- 3 or 4 state model:
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
mod <- two_facilitation_model()
mod

# Show parameters:
parms(mod)

# Tweak parameters:
parms(mod)["g"] <- 0.2
parms(mod)["gamma1"] <- 0.10
parms(mod)["u"] <- 5 
times(mod) <- c(from = 0, to = 1000, by = 1)
```

- Launch a simulation and plot the result:

```
mod_run <- sim(mod)
plotnp(mod_run)
```

## Repo organization

Here is the tree of the directory:

```
.
├── draft
│   ├── analytic
│   ├── code
│   │   └── code_hp_windows
│   ├── results
│   └── simulation
│       ├── CA_mod3
│       │   └── CA_mod3
│       ├── CA.mod4
│       ├── model1
│       ├── model2
│       └── pairs_mod3
├── inst
│   ├── figs
│   │   └── four_states
│   └── r
├── man
├── R
├── tests
│   └── testthat
└── vignettes
```

The repo is organized as a R package. The functions are defined in the `R`
folder and the results are reported in the `vignettes`.

### Functions

```
.
├── analysis_methods.R 
├── ca_4_states_model.R
├── ca_transitions_functions.R
├── cellular_automata_template.R
├── check_arguments.R
├── deprecated.R
├── indirect_facilitation.R
├── misc.R
├── model_template.R
├── ode_template.R
├── pa_3_states_model.R
├── pa_4_states_model.R
├── pa_transition_functions.R
├── plot_methods.R
└── running_methods.R
```

