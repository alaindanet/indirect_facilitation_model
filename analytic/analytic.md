---
title: 'Thesis, Chap 1: Modelisation'
output: html_document
---

# Objectives:

1. Find more coexistence between the two species.
2. Get more spatial autocorrelation between the 2 species.

# Tasks:

* Find analytical solution by pair approximation
* Find a smart way to analyse model.

# First step: simple model with 3 states

* We have 3 states ($+$, $0$, $-$).
* We have six different pairs ($\rho_{++}$, $\rho_{+0}$, $\rho_{+-}$, $\rho_{00}$, $\rho_{0-}$, $\rho_{--}$) since $\rho_{\sigma\sigma'} = \rho_{\sigma'\sigma}$.
* We have four conservation equations:
  1. $\rho_{+} + \rho_{0} + \rho_{-} = 1$
  2. $\rho_{+} = \rho_{++} + \rho_{+0} + \rho_{+-}$
  3. $\rho_{0} = \rho_{00} + \rho_{+0} + \rho_{0-}$
  4. $\rho_{-} = \rho_{--} + \rho_{-0} + \rho_{+-}$
* There are 3 singleton variables: $\rho_{+}$, $\rho_{0}$, $\rho_{-}$

So we need $6+3-4 = 5$ equations to solve this system. As in Kéfi et al. (2007), I choose $\frac{d\rho_{++}}{dt}$, $\frac{d\rho_{+-}}{dt}$, $\frac{d\rho_{--}}{dt}$, $\frac{d\rho_{+}}{dt}$, $\frac{d\rho_{-}}{dt}$.

## 1. $\frac{d\rho_{++}}{dt}$

$\frac{d\rho_{++}}{dt} = 2\rho_{0+}w_{0,+} - 2\rho_{++}w_{+,0}$

With $\rho_{+0} = \rho_+ - \rho_{++} - \rho_{0-}$, $q_{+|0} = \frac{\rho_{+0}}{\rho_0} = \frac{\rho_+ - \rho_{++} - \rho_{0-}}{\rho_0}$.

We have:

$$
\begin{aligned}
\frac{d\rho_{++}}{dt} =& 2\rho_{+0} [ (\delta\rho_{+} + \frac{z-1}{z}(1-\delta)\frac{\rho_{+0}}{\rho_0}) \times (b - c\rho_+ - g(1- \frac{z-1}{z}*\frac{\rho_{+0}}{\rho_{0} } n)) ] - 2\rho_{++}m \\
\frac{d\rho_{++}}{dt} =& 2(\rho_+ - \rho_{++} - \rho_{0-})[(\delta\rho_{+} + \frac{z-1}{z}(1-\delta)\frac{\rho_+ - \rho_{++} - \rho_{+-}}{1-\rho_- - \rho_+}) \\
& \times(b - c\rho_+ - g(1- \frac{z-1}{z}*\frac{\rho_+ - \rho_{++} - \rho_{+-}}{1-\rho_- - \rho_+}n))] \\
& - 2\rho_{++}m
\end{aligned}
$$

Solve: (je n'ai pas réussi à simplifier des termes après développement pour celle ci)

$$
\begin{aligned}
\frac{d\rho_{++}}{dt} =& 0 \\
2\rho_{+0} [ (\delta\rho_{+} + \frac{z-1}{z}(1-\delta)\frac{\rho_{+0}}{\rho_0}) \times (b - c\rho_+ - g(1- \frac{z-1}{z}*\frac{\rho_{+0}}{\rho_{0} } n)) ] - 2\rho_{++}m =& 0 \\
 2\rho_{+0} [ (\delta\rho_{+} + \frac{z-1}{z}(1-\delta)\frac{\rho_{+0}}{\rho_0}) \times (b - c\rho_+ - g(1- \frac{z-1}{z}*\frac{\rho_{+0}}{\rho_{0} } n)) ] =& 2\rho_{++}m \\
 2\rho_{+0} [ \delta\rho_{+}b - \delta\rho_{+}^2c - \delta\rho_{+}g(1- \frac{z-1}{z}*\frac{\rho_{+0}}{\rho_{0} } n) + b\frac{z-1}{z}(1-\delta)\frac{\rho_{+0}}{\rho_0} - c\rho_+\frac{z-1}{z}(1-\delta)\frac{\rho_{+0}}{\rho_0} - g(1- \frac{z-1}{z}*\frac{\rho_{+0}}{\rho_{0} } n)\frac{z-1}{z}(1-\delta)\frac{\rho_{+0}}{\rho_0} ] =& 2\rho_{++}m \\
 2\rho_{+0} [ \delta\rho_{+}b - \delta\rho_{+}^2c - \delta\rho_{+}g + \delta\rho_{+}g\frac{z-1}{z}*\frac{\rho_{+0}}{\rho_{0} }n - c\rho_+\frac{z-1}{z} + c\delta\rho_+\frac{z-1}{z} + b\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0} - b\delta\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0} - (g - g\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0}n)(\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0}(1-\delta)) ] =& 2\rho_{++}m \\
 \rho_{+0} [ \delta\rho_{+}b - \delta\rho_{+}^2c - \delta\rho_{+}g + \delta\rho_{+}g\frac{z-1}{z}*\frac{\rho_{+0}}{\rho_{0} }n - c\rho_+\frac{z-1}{z} + c\delta\rho_+\frac{z-1}{z} + b\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0} - b\delta\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0} - (g - g\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0}n)(\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0}-\delta\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0} ) ] =& \rho_{++}m \\
 \rho_{+0} [ \delta\rho_{+}b - \delta\rho_{+}^2c - \delta\rho_{+}g + \delta\rho_{+}g\frac{z-1}{z}*\frac{\rho_{+0}}{\rho_{0} }n - c\rho_+\frac{z-1}{z} + c\delta\rho_+\frac{z-1}{z} + b\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0} - b\delta\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0} - (g\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0} - g\delta\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0} - g(\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0})^2 + g\delta(\frac{z-1}{z}\frac{\rho_{+0}}{\rho_0})^2) =& \rho_{++}m \\
\frac{ \rho_{+0}(\delta\rho_{+} + \frac{z-1}{z}(1-\delta)\frac{\rho_{+0}}{\rho_0}) \times (b - c\rho_+ - g(1- \frac{z-1}{z}*\frac{\rho_{+0}}{\rho_{0} } n))}{m} =& \rho_{++} 
\end{aligned}
$$

## 2. $\frac{d\rho_{+-}}{dt}$

$\frac{d\rho_{+-}}{dt} = \rho_{0+}w_{0,-} + \rho_{0-}w_{0,+} - \rho_{+-}(w_{-,0}+w_{+,0})$

With: $\rho_{0-} = \rho_- - \rho_{--} - \rho_{+-}$, $q_{+|-} = \frac{\rho_{+-}}{\rho_-}$.

We have: 

$$
  \begin{aligned}
  \frac{d\rho_{+-}}{dt} =& d(\rho_+ - \rho_{++} - \rho_{+-}) \\
  &+ (\rho_- - \rho_{--} - \rho_{-+})[(\delta\rho_{+} + \frac{z-1}{z}(1-\delta)\frac{\rho_+ - \rho_{++} - \rho_{+-}}{1-\rho_--\rho_+}) \\
  & \times(b - c\rho_+ - g(1- \frac{z-1}{z}*\frac{\rho_+ - \rho_{++} - \rho_{+-}}{1-\rho_- - \rho_+}n))] \\
  &- \rho_{+-}(m + r + \frac{f}{z} + \frac{z-1}{z}\frac{\rho_{+-}}{\rho_-}f)
  \end{aligned}
$$

Solve: 


$$
\begin{aligned}
\frac{d\rho_{+-}}{dt} =& 0 \\
- \rho_{+-}(m + r + \frac{f}{z} + \frac{z-1}{z}\frac{\rho_{+-}}{\rho_-}f) =& - d\rho_{+0} - \rho_{0-}[(\delta\rho_{+} + \frac{z-1}{z}(1-\delta)\frac{\rho_{+0}}{\rho_0}) \times(b - c\rho_+ - g(1- \frac{z-1}{z}*\frac{\rho_{+0}}{\rho_0}n))] \\
\rho_{+-} =&  \frac{ d\rho_{+0} + \rho_{0-}[(\delta\rho_{+} + \frac{z-1}{z}(1-\delta)\frac{\rho_{+0}}{\rho_0}) \times(b - c\rho_+ - g(1- \frac{z-1}{z}*\frac{\rho_{+0}}{\rho_0}n))] }{m + r + \frac{f}{z} + \frac{z-1}{z}\frac{\rho_{+-}}{\rho_-}f} \\
\end{aligned}
$$

## 3. $\frac{d\rho_{--}}{dt}$

$\frac{d\rho_{--}}{dt} = 2\rho_{0-}w_{0,-} - 2\rho_{--}w_{-,0}$

With: $\rho_{0-}=\rho_{-} - \rho_{--} - \rho_{+-}$, $w_{0,-}=d, q_{+|-} = \frac{\rho_{+-}}{\rho_-}$, $w_{-,0}= r + \frac{z-1}{z}\frac{\rho_{+-}}{\rho_{-}}f$

We have:

$$
\frac{d\rho_{--}}{dt} = 2(\rho_- - \rho_{--} - \rho_{+-})d - 2\rho_{--}(r + \frac{z-1}{z}\frac{\rho_{+-}}{\rho_-}f)
$$

Solve:

$$
\begin{aligned}
\frac{d\rho_{--}}{dt} =& 0 \\
2(\rho_- - \rho_{--} - \rho_{+-})d - 2\rho_{--}(r + \frac{z-1}{z}\frac{\rho_{+-}}{\rho_-}f) =& 0 \\
2d\rho_- - 2d\rho_{--} - 2d\rho_{+-} - 2r\rho_{--} + 2\frac{z-1}{z}\frac{\rho_{+-}}{\rho_-}f\rho_{--} =& 0 \\
-2\rho_{--}(d + r + \frac{z-1}{z}\frac{\rho_{+-}}{\rho_-}f) =& 2d\rho_{+-} - 2d\rho_- \\
- \rho_{--}(d + r + \frac{z-1}{z}\frac{\rho_{+-}}{\rho_-}f) =& d(\rho_{+-} - \rho_-) \\
=& d(-\rho_{--} - \rho_{-0}) \\
d\rho_{--} - \rho_{--}(d + r + \frac{z-1}{z}\frac{\rho_{+-}}{\rho_-}f) =& - d\rho_{-0} \\
- \rho_{--}(r + \frac{z-1}{z}\frac{\rho_{+-}}{\rho_-}f) =& - d\rho_{-0} \\
\rho_{--} =& \frac{ d\rho_{-0} }{ r + \frac{z-1}{z}\frac{\rho_{+-}}{\rho_-}f }
\end{aligned}
$$

## 4. $\frac{d\rho_{+}}{dt}$

$\frac{d\rho_{+}}{dt}= \rho_{0}w_{0,+}$

With: $\rho_0 = 1- \rho_+ - \rho_- $, $w_{0,+} - \rho_+m $

We have:

$$
\begin{aligned}
\frac{d\rho_{+}}{dt}=& (1- \rho_+ - \rho_-)[(\delta\rho_{+} + (1-\delta)\frac{\rho_+ - \rho_{++} - \rho_{+-}}{1-\rho_- - \rho_+}) \\
  &\times(b - c\rho_+ - g(1- \frac{\rho_+ - \rho_{++} - \rho_{+-}}{1-\rho_- - \rho_+}n))] \\
  & - \rho_+m
\end{aligned} 
$$

Solve:

$$
\begin{aligned}
\frac{d\rho_{+}}{dt} =& 0 \\
\rho_0(\delta\rho_{+} + (1-\delta)\frac{\rho_{+0}}{\rho_0}) \times(b - c\rho_+ - g(1- \frac{\rho_{+0}}{\rho_0}n))- \rho_+m =& 0 \\
\frac{ \rho_0(\delta\rho_{+} + (1-\delta)\frac{\rho_{+0}}{\rho_0}) \times(b - c\rho_+ - g(1- \frac{\rho_{+0}}{\rho_0}n)) }{m} =& \rho_+
\end{aligned} 
$$

## 5. $\frac{d\rho_{-}}{dt}$

$\frac{d\rho_{-}}{dt} = \rho_{0}w_{0,-} - \rho_{-}w_{-,0}$

We have:

$$
\begin{aligned}
\frac{d\rho_{-}}{dt} = d(1-\rho_{+} - \rho_{-}) - (r + \frac{\rho_{+-}}{\rho_-}f)\rho_{-}
\end{aligned}
$$

Solve:

$$
\begin{aligned}
\frac{d\rho_{-}}{dt} =& 0 \\
d(1-\rho_{+} - \rho_{-}) - (r + \frac{\rho_{+-}}{\rho_-}f)\rho_{-} =& 0 \\
d - d\rho_{+} - d\rho_{-} - r\rho_{-} - \frac{\rho_{+-}}{\rho_-}f\rho_{-} =& 0 \\
d - d\rho_{+} - \rho_{-}(d + r) - \rho_{+-}f =& 0 \\
d\rho_{0} + d\rho_{-} - \rho_{-}(d + r) =& \rho_{+-}f \\
-\rho_{-}r =& \rho_{+-}f - d\rho_{0} \\
\rho_{-} =& \frac{ -\rho_{+-}f + d\rho_{0} }{r} \\
\end{aligned}
$$
