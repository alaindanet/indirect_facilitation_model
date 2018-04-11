# Effect of indirect and direct facilitation on coexistence and stability of a dryland ecosystem  

## TODO  

- [x] Ecological model
- [x] ODE system function 
- [x] ODEs  
- [x] Transition equation functions  
- [x] Test functions
- [/] Document functions
- [ ] Framework to build simulation diagrams
    - [x] Run simulations over a range of parameters  
    - [x] Extract the average of the last timesteps
    - [ ] Check if Transition probabilities are negative (return a WARNING ?)
    - [ ] Ajouter un seuil d'extinction
    - [ ] Check if simulations have reached stability
    - [ ] Define a plot method 

## Objectives

Develop a theoretical model of two species in competition and add an interaction
of indirect facilitation. A nurse species is protected against grazing can protect the saplings of a protégée species by associative protection. We make the simple assumption of a growth-defence trade-off.

The cells can take 3 states: occupied by a nurse, a protégée or be empty.

## Three states model - The nurse, the empty & the protégée 


* We have 3 states ($N$, $P$, $0$).
* We have six differeNt pairs ($\rho_{NN}$, $\rho_{N0}$, $\rho_{NP}$, $\rho_{PP}$, $\rho_{P0}$, $\rho_{00}$) since $\rho_{\sigma\sigma'} = \rho_{\sigma'\sigma}$.

* We have four coNservatioN equatioNs:
  1. $\rho_{N} + \rho_{0} + \rho_{P} = 1$
  2. $\rho_{N} = \rho_{NN} + \rho_{N0} + \rho_{NP}$
  3. $\rho_{0} = \rho_{00} + \rho_{N0} + \rho_{P0}$
  4. $\rho_{P} = \rho_{PP} + \rho_{P0} + \rho_{NP}$
* There are 3 singleton variables: $\rho_{N}$, $\rho_{0}$, $\rho_{P}$

So we need $6+3-4 = 5$ equations to solve this system:

 * $\frac{d\rho_{NN}}{dt}$
 * $\frac{d\rho_{NP}}{dt}$
 * $\frac{d\rho_{PP}}{dt}$
 * $\frac{d\rho_{N}}{dt}$
 * $\frac{d\rho_{P}}{dt}$

### Transition equations

$$
w_{0,N|P} = (\delta \rho_{N} + (1-\delta) \frac{(z-1)}{z} \times \frac{ \rho_{N0}
}{\rho_{0}}) \times (b - c\rho_{+} - \gamma)
$$

$$
w_{0,N|N} = (\delta \rho_{N} + \frac{(1-\delta)}{z} +  (1-\delta) \frac{(z-1)}{z} \frac{\rho_{N0}}{\rho_{0}}) \times (b - c\rho_{+} - \gamma)
$$

$$
w_{0,N} = (\delta \rho_{N} + (1-\delta) \frac{\rho_{N0}}{\rho_{0}}) \times (b - c\rho_{+} - \gamma)
$$

$$
w_{0,P} = (\delta \rho_{P} + (1 - \delta) \frac{\rho_{P0}}{\rho_{0}})
\times (b - c\rho_{+} - g(1 - \frac{\rho_{N0}}{\rho_{0}}f(\gamma)))
$$

$$
w_{0,P|N} = (\delta\rho_{P} + (1 - \delta) \frac{(z-1)}{z}
\frac{\rho_{P0}}{\rho_{0}}) \times (b - c\rho_{+} - g(1 -
(\frac{(z-1)}{z} \frac{\rho_{N0}}{\rho_{0}} + \frac{1}{z})f(\gamma)))
$$

$$
w_{0,P|P} = (\delta\rho_{P} + (1 - \delta) \frac{(z-1)}{z}
\frac{\rho_{P0}}{\rho_{0}}) \times (b - c\rho_{+} - g(1 -
\frac{(z-1)}{z} \frac{\rho_{N0}}{\rho_{0}}f(\gamma)))
$$

$$
w_{0,P} = (\delta\rho_{P} + (1 - \delta) \frac{\rho_{P0}}{\rho_{0}}) \times (b -
c\rho_{+} - g(1 - \frac{\rho_{N0}}{\rho_{0}}f(\gamma)))
$$

$$
w_{P,0} = w_{P,0} = m
$$

### ODEs system 

WARNING: I am not sure of the 2 times

$$
\begin{aligned}
\frac{d\rho_{NP}}{dt} = \rho_{P0} \times w_{0,N|P} + \rho_{N0} \times w_{0,P|N}
- \rho_{NP} \times (w_{P,0} + w_{N,0})
\end{aligned}
$$

$$
\begin{aligned}
\frac{d\rho_{PP}}{dt} = 2\rho_{P0} \times w_{0,P|P} - 2\rho_{PP} \times w_{P,0}\\
\end{aligned}
$$

$$
\begin{aligned}
\frac{d\rho_{NN}}{dt}= 2\rho_{N0} \times w_{0,N|N} - 2\rho_{NN} \times w_{N,0}  
\end{aligned}
$$

$$
\begin{aligned}
\frac{d\rho_{P}}{dt} = \rho_{0} \times w_{0,P} - \rho_{P} \times w_{P,0}
\end{aligned}
$$

$$
\begin{aligned}
\frac{d\rho_{N}}{dt} = \rho_{0} \times w_{0,N} - \rho_{N} \times w_{N,0}
\end{aligned}
$$

### Definition of the associative protection $n = f(\gamma)$  

$n(\gamma) = $ 

Le Galliard et al. (2003) models the cost of altruism as: 

$$ 
\gamma = Kn^{u}
$$

with $n$, the level of altruism, $u$, the shape of the relationship between the
associative protection and the cost.

Lion et al. (2007) used the same kind of relationship but assume that the cost
was dependant of the environment (i.e. of the neighboring cells).

 
In our case, we can assume that there is a relationship between the associative
protection and the cost of protection.

Let $n(\gamma) = kn^{u}$
