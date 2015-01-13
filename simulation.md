---
title: "Simulation"
output: html_document
---

## Tasks

1. Find more coexistence
2. Find more spatial correlation between species


## Transitions rates

### Colonisation
$$
  w_{ \left\{0,+_N \right\} } = \left( \delta_N\rho_{+_N} + \left( 1 - \delta_N \right)q_{+_N|0}\right) \left(b_N-c_N\rho_{+_N} -  c_{PN}\rho_{+_P} - gq_{+_P|+_N}p \right)
$$
$$
  w_{ \left\{0,+_P \right\} } = \left( \delta_P\rho_{+_P} + \left( 1 - \delta_P \right)q_{+_P|0}\right) \left(b_P-c_P\rho_{+_P} -  c_{NP}\rho_{+_N} - g(1 - q_{+_N|+_P}n) \right)
$$

### Regeneration

$$
w_{ \left\{-,0 \right\} } = r + q_{+|-} f
$$

### Mortality

$$ w_{ \left\{ +_N,0 \right\} }  = m_N $$ $$w_{ \left\{ +_P,0 \right\} }  = m_P$$

### Degradation

$$w_{ \left\{0,- \right\} } = d$$

## Basic parameter values

| Parameters  | Values  |
| ----------- | ------- |
| $        |
|            |        |
|            |        |

## 1. Coexistence

To get more coexistence between our 2 species, we have to explore  inter and intra competition parameters. There is 2 main scenarii: 1) Intra competition is high and inter is low. 2) Intra competition is low and inter is high as in my previous work.

I defined presence of a species when its cover is higher than 5%.

### 1.1 Scenario 1 : High intra competition and low inter

##### Parameters

$c_1 = c_2 = 0.4$, $c_{12} = [0.1,0.2,0.3,0.4]$, $c_{21} = [0.1,0.2,0.3,0.4]$, $b=[0.5,0.6,0.7]$, $g=[0,0.05,0.1]$, $d=0.1$, $r=0.0$, $m=0.2$, $\delta = 0.1$, $n =[0,1]$

##### Results

#### Scenario 1.1: Intra competition higher than inter: Competition lower than previously

##### Parameters

$c_1 = c_2 = 0.2$, $c_{12} = [0.05,0.1,0.15,0.2]$, $c_{21} = [0.05,0.1,0.15,0.2]$, $b=[0.5,0.6,0.7]$, $g=[0,0.05,0.1]$, $d=0.1$, $r=0.0$, $m=0.2$, $\delta = 0.1$, $n =[0,1]$

##### Results

Coexistence occured but the clustering between the 2 species were always lower than 1. The nursing has not been very effective (26 when n=0 vs 31 when n=1). My point of view was two species can potentially coexist is there are more intra- than inter-competition. But stay, these species can so coexist whithout nursing effect and the clustering is lower than 1. So which hypothesis can make the cooccurence between my two species higher ?

### 1.2 Global competition and local intercompetition


