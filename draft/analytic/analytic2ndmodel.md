---
title: 'Thesis, Chap 1: Modelisation: 3 states model'
output: html_document
---

# Objectives:

1. Find more coexistence between the two species.
2. Get more spatial autocorrelation between the 2 species.

# Tasks:

* Find analytical solution by pair approximation
* Find a smart way to analyse model.

# First step: simple model with 3 states

* We have 3 states ($n$, $p$, $0$).
* We have six different pairs ($\rho_{nn}$, $\rho_{n0}$, $\rho_{np}$, $\rho_{pp}$, $\rho_{0p}$, $\rho_{00}$) since $\rho_{\sigma\sigma'} = \rho_{\sigma'\sigma}$.

* We have four conservation equations:
  1. $\rho_{n} + \rho_{0} + \rho_{p} = 1$
  2. $\rho_{n} = \rho_{nn} + \rho_{n0} + \rho_{np}$
  3. $\rho_{0} = \rho_{00} + \rho_{n0} + \rho_{0p}$
  4. $\rho_{p} = \rho_{pp} + \rho_{p0} + \rho_{np}$
* There are 3 singleton variables: $\rho_{n}$, $\rho_{0}$, $\rho_{p}$

So we need $6+3-4 = 5$ equations to solve this system:
* $\frac{d\rho_{nn}}{dt}$
* $\frac{d\rho_{np}}{dt}$
* $\frac{d\rho_{pp}}{dt}$
* $\frac{d\rho_{n}}{dt}$
* $\frac{d\rho_{p}}{dt}$

Equations:

