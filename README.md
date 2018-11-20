
<!-- README.md is generated from README.Rmd. Please edit that file -->

# indirectfacilitation

The goal of indirectfacilitation is to …

## Installation

You can install indirectfacilitation from github with:

``` r
# install.packages("devtools")
devtools::install_github("alaindanet/indirect_facilitation_model")
```

## Example

This is a basic example which shows you how to solve a common problem:

  - Load the packages:

<!-- end list -->

    devtools::load_all()
    
    library(simecol)
    library(ggplot2)
    library(tibble)
    library(magrittr)

  - Load the model and specify the parameters:

<!-- end list -->

    mod <- two_facilitation_model()
    mod
    
    # Show parameters:
    parms(mod)
    
    # Tweak parameters:
    parms(mod)["g"] <- 0.2
    parms(mod)["gamma1"] <- 0.10
    parms(mod)["u"] <- 5 
    times(mod) <- c(from = 0, to = 1000, by = 1)

  - Launch a simulation and plot the result:

<!-- end list -->

    mod_run <- sim(mod)
    plotnp(mod_run)

## Todo

### Paper

  - \[ \] Fig 1: methods
      - \[ \] Add temporal dynamic CA + PA
  - \[ \] Fig 2: multi-states
      - \[ \] Add arrows
  - \[x\] Fig 3: bifurcation details
  - \[x\] Fig 4: Clustering cellular automata, C++ et Cnp
    \(f(\gamma, u)\) and \(f(\gamma, u)\)
      - \[ \] (Change space between
        panel)\[<https://github.com/tidyverse/ggplot2/issues/2464#issuecomment-384980197>\]
      - \[ \] (Change margin of the plot for
        plot\_grid)\[<https://github.com/wilkelab/cowplot/issues/31#issuecomment-292748114>\]
      - \[ \] Fix B and E legend decrease a bit the lowest value (some
        points disapeared in fig B right panel, en haut à gauche)

### Run

1.  \[ \] Run pair-approximation model with automated end of run
2.  \[ \] Replace the file with `use_data()`

### Implementation

  - \[ \] Clustering:
      - \[ \] plot NP
      - \[ \] plot (x2 - x1) / x1 for comparison CA/PA
      - \[ \] run CA for gradient facilitation/dispersion
  - \[ \] plot: separate lines (create groups: e.g: together & \<
    threshold, together & \>= threshold, low\_together …)
      - \[ \] Create groups
      - \[ \] Plot groups (linetype problem)
  - \[x\] Check if simulations have reached stability
      - \[x\] Implementation of a custom solver (“lsodar”)
      - \[x\] Test the implementation
  - \[ \] sim\_multi():
      - \[ \] Save a simecol object
      - \[ \] Save a time argument
      - \[ \] Save init values
      - \[ \] Save baseline parameters
      - \[ \] Save param\_combination
      - \[ \] param\_combination: matrix of parameter combination
      - \[ \] be able to sim from a sim\_multi object
  - \[ \] Generalize run\_scenarii\_gradient and run\_2d\_gradient:
    f(b,g)
      - \[ \] check\_param(model\_spec, names(gradient)))
      - \[ \] Add check for argument conformity (check\_param,
        check\_inits) of valid starting values (scenario = “custom”), or
        a named list (provided by init\_scenarii())
      - \[ \] Add methods: init.scenarii(), param.scenarii() (S4
        methods)
        <https://stackoverflow.com/questions/12100856/combining-s4-and-s3-methods-in-a-single-function>
      - \[ \] See [here](https://stackoverflow.com/a/23863893/5968131)
  - \[x\] Compare the effect of dispersal, facilitation strength,
    paturâge, aridité on clustering  
  - \[ \] Document functions
