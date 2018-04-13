# Effect of indirect and direct facilitation on coexistence and stability of a dryland ecosystem  

## TODO  

- [x] Ecological model
- [x] ODE system function 
- [x] ODEs  
- [x] Transition equation functions  
- [x] Test functions
- [/] Document functions
- [/] Framework to build simulation diagrams
    - [x] Run simulations over a range of parameters  
    - [x] Extract the average of the last timesteps
    - [x] Check if simulations become weird (return a WARNING ?)
    - [x] Ajouter un seuil d'extinction
    - [ ] Check if simulations have reached stability
    - [x] Define a plot method 
    - [ ] Generalize functions for any parameters 

## Objectives

Develop a theoretical model of two species in competition and add an interaction
of indirect facilitation. A nurse species is protected against grazing can protect the saplings of a protégée species by associative protection. We make the simple assumption of a growth-defence trade-off.

The cells can take 3 states: occupied by a nurse, a protégée or be empty.

## Workflow 

```
library(indirectfacilitation)

```
