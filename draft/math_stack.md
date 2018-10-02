# Conditional probabilities and pair approximation model

Let's consider a square lattice of cells. Each cell can be either occupied by a species (1 or 2) or be empty (0).
Each cell can be either in state 1, 2 or 0.  

In the pair approximation model, I would like to compute the clustering of the
species, i.e. the clustering of the occupied cells ($+$). It is define as:

* $C_{++} = \frac{q_{+|+}}{\rho_{+}} = \frac{\rho_{++}}{\rho_{+}^2}$

Where

* $q_{+|+}$ is the conditional probability to find an occupied cell in the surrounding cells knowing that the focal cell is occupied
* $\rho_+$ is the density of occupied cells in the landscape, defines as: $\rho_+ = \rho_1 + \rho_2$ 
* $\rho_++$ is the density of occupied cell **pairs** in the landscape, defines as: $\rho_{++} = \rho_{11} + \rho_{12} + \rho_{21} + \rho_{22}$ 

* $q_{i|j} = \frac{\rho_{ij}}{\rho_{i}}$

One approach successfully describes the clustering but the second does not and I
do not find why.

## A first approach (Works)

* $q_{+|+} = \frac{\rho_{++}}{\rho_{+}}$

* $\rho_{++} = \rho_{11} + \rho_{12} + \rho_{21} + \rho_{22}$

Knowing that $\rho_{12} = \rho_{21}$:

* $q_{+|+} = \frac{\rho_{11} + 2\rho_{12} + \rho_{22}}{\rho_{1} + \rho_{2}}$

Hence:

$$
C_{++} = \frac{\rho_{11} + 2\rho_{12} + \rho_{22}}{\rho_{1} + \rho_{2}} \times \frac{1}{\rho_{+}} \\
C_{++} = \frac{\rho_{11} + 2\rho_{12} + \rho_{22}}{(\rho_{1} + \rho_{2})^2}
$$

This approach works, I have checked it by running simulation of cellular
automata.

## A second approach (Do not works)

Let's define $q_{+|+}$.

It is the probability for one cell of state $1$ to be surrounded by occupied
cells +  the probability for one cell of state $2$ to be surrounded by occupied
cells. It means that:

$$
q_{+|+} = q_{+|1} + q_{+|2}
$$

and $q_{+|1}$ and $q_{+|2}$ can be defined as:

* $q_{+|1} = q_{1|1} + q_{2|1}$
* $q_{+|2} = q_{1|2} + q_{2|2}$

So: 

\begin{align}
q_{+|+} & = q_{1|1} + q_{2|1} + q_{1|2} + q_{2|2} \\
 & = \frac{\rho_{11}}{\rho_{1}} \frac{\rho_{12}}{\rho_{1}} + \frac{\rho_{21}}{\rho_{2}} + \frac{\rho_{22}}{\rho_{2}} \\
 & = \frac{\rho_{11} + \rho_{12}}{\rho_{1}} + \frac{\rho_{21} + \rho_{22}}{\rho_{2}} \\ 
 & = \frac{\rho_{2}(\rho_{11} + \rho_{12})}{\rho_{1}\rho_{2}} + \frac{\rho_{1}(\rho_{21} + \rho_{22})}{\rho_{1}\rho_{2}} \\ 
 & = \frac{\rho_{2}(\rho_{11} + \rho_{12}) + \rho_{1}(\rho_{21} + \rho_{22}) }{\rho_{1}\rho_{2}} \\
\end{align}

Obviously this formulation of $q_{+|+}$ is different for the first one. In
simulation, this formulation gives $C_{++}$ twice higher than with a first
approach.

I do not what is the mistake in this approach but I guess that the following
assertion is false: 

$$q_{+|+} & = q_{1|1} + q_{2|1} + q_{1|2} + q_{2|2}$$

I would be very helpful for me to know why.

