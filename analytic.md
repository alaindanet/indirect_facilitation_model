---
title: 'Thesis, Chap 1: Modelisation'
output: html_document
---

## Objectives:

1. Find more coexistence between the two species.
2. Get more spatial autocorrelation between the 2 species.

## Tasks:

* Find analytical solution by pair approximation
* Find a smart way to analyse model.

## I. Simple model with 3 states

* We have 3 states ($+$, $0$, $-$).
* We have six different pairs ($\rho_{++}$, $\rho_{+0}$, $\rho_{+-}$, $\rho_{00}$, $\rho_{0-}$, $\rho_{--}$) since $\rho_{\sigma\sigma'} = \rho_{\sigma'\sigma}$.
* We have four conservation equations:
  1. $\rho_{+} + \rho_{0} + \rho_{-} = 1$
  2. $\rho_{+} = \rho_{++} + \rho_{+0} + \rho_{+-}$
  3. $\rho_{0} = \rho_{00} + \rho_{+0} + \rho_{0-}$
  4. $\rho_{-} = \rho_{--} + \rho_{-0} + \rho_{+-}$
* There are 3 singleton variables: $\rho_{+}$, $\rho_{0}$, $\rho_{-}$

So we need $6+3-4 = 5$ equations to solve this system. As in Kéfi et al. (2007), I choose $\frac{d\rho_{++}}{dt}$, $\frac{d\rho_{+-}}{dt}$, $\frac{d\rho_{--}}{dt}$, $\frac{d\rho_{+}}{dt}$, $\frac{d\rho_{-}}{dt}$.

### Write equations

#### 1. $\frac{d\rho_{++}}{dt}$

$\frac{d\rho_{++}}{dt} = 2\rho_{0+}w_{0,+} - 2\rho_{++}w_{+,0}$

With $\rho_{+0} = \rho_+ - \rho_{++} - \rho_{0-}$, $q_{+|0} = \frac{\rho_{+0}}{\rho_0} = \frac{\rho_+ - \rho_{++} - \rho_{0-}}{\rho_0}$.

We have:

$$
\begin{aligned}
\frac{d\rho_{++}}{dt} =& 2(\rho_+ - \rho_{++} - \rho_{0-})[(\delta\rho_{+} + \frac{z-1}{z}(1-\delta)\frac{\rho_+ - \rho_{++} - \rho_{+-}}{1-\rho_- - \rho_+}) \\
& \times(b - c\rho_+ - g(1- \frac{z-1}{z}\times\frac{\rho_+ - \rho_{++} - \rho_{+-}}{1-\rho_- - \rho_+}n))] \\
& - 2\rho_{++}m
\end{aligned}
$$


#### 2. $\frac{d\rho_{+-}}{dt}$

$\frac{d\rho_{+-}}{dt} = \rho_{0+}w_{0,-} + \rho_{0-}w_{0,+} - \rho_{+-}(w_{-,0}+w_{+,0})$

With: $\rho_{0-} = \rho_- - \rho_{--} - \rho_{+-}$, $q_{+|-} = \frac{\rho_{+-}}{\rho_-}$.

We have: 

$$
  \begin{aligned}
  \frac{d\rho_{+-}}{dt} =& d(\rho_+ - \rho_{++} - \rho_{+-}) \\
  &+ (\rho_- - \rho_{--} - \rho_{-+})[(\delta\rho_{+} + \frac{z-1}{z}(1-\delta)\frac{\rho_+ - \rho_{++} - \rho_{+-}}{1-\rho_--\rho_+}) \\
  & \times(b - c\rho_+ - g(1- \frac{z-1}{z}\times\frac{\rho_+ - \rho_{++} - \rho_{+-}}{1-\rho_- - \rho_+}n))] \\
  &- \rho_{+-}(m + r + \frac{f}{z} + \frac{z-1}{z}\frac{\rho_{+-}}{\rho_-}f)
  \end{aligned}
$$

#### 3. $\frac{d\rho_{--}}{dt}$

$\frac{d\rho_{--}}{dt} = 2\rho_{0-}w_{0,-} - 2\rho_{--}w_{-,0}$

With: $\rho_{0-}=\rho_{-} - \rho_{--} - \rho_{+-}$, $w_{0,-}=d, q_{+|-} = \frac{\rho_{+-}}{\rho_-}$, $w_{-,0}= r + \frac{z-1}{z}\frac{\rho_{+-}}{\rho_{-}}f$

We have:

$$
\frac{d\rho_{--}}{dt} = 2(\rho_- - \rho_{--} - \rho_{+-})d - 2\rho_{--}(r + \frac{z-1}{z}\frac{\rho_{+-}}{\rho_-}f)
$$


#### 4. $\frac{d\rho_{+}}{dt}$

$\frac{d\rho_{+}}{dt}= \rho_{0}w_{0,+} $

With: $\rho_0 = 1- \rho_+ - \rho_- $, $w_{0,+} - \rho_+m$

We have:

$$
\begin{aligned}
\frac{d\rho_{+}}{dt}=& (1- \rho_+ - \rho_-)[(\delta\rho_{+} + (1-\delta)\frac{\rho_+ - \rho_{++} - \rho_{+-}}{1-\rho_- - \rho_+}) \\
  &\times(b - c\rho_+ - g(1- \frac{\rho_+ - \rho_{++} - \rho_{+-}}{1-\rho_- - \rho_+}n))] \\
  & - \rho_+m
\end{aligned}
$$

#### 5. $\frac{d\rho_{-}}{dt}$

$\frac{d\rho_{-}}{dt} = \rho_{0}w_{0,-} - \rho_{-}w_{-,0}$

We have:

$$
\begin{aligned}
\frac{d\rho_{-}}{dt} = d(1-\rho_{+} - \rho_{-}) - (r + \frac{\rho_{+-}}{\rho_-}f)\rho_{-}
\end{aligned}
$$

### Solve Equations for $\frac{dp_{\sigma}}{dt} = 0$ and $\frac{dp_{\sigma\sigma}}{dt} = 0$

#### 1. $\frac{d\rho_{++}}{dt} = 0$

#### 2. $\frac{d\rho_{+-}}{dt} = 0$

#### 3. $\frac{d\rho_{--}}{dt} = 0$

#### 4. $\frac{d\rho_{+}}{dt} = 0$

#### 5. $\frac{d\rho_{-}}{dt} = 0$

$$
\begin{aligned}
d(1-\rho_+-\rho_-) - (r+\frac{\rho_{+-}}{\rho_-}f)\rho_- =& 0 \\
d(1- \rho_+) - d\rho_- - r\rho_- - f\rho_-\frac{\rho_{+-}}{\rho_-} =& 0 \\
d\rho_0 + d\rho_- - d\rho_- - r\rho_- - f\rho_{+-} =& 0 \\
- r\rho_- =& -d\rho_0 + f\rho_{+-} \\
r\rho_- =& d\rho_0 - f\rho_{+-} \\
\rho_- =& \frac{ d\rho_0 - f\rho_{+-}}{r} \\
\end{aligned}
$$