#############################
#  Test for the clustering  #
#############################

# Variables 
## Rho: n,p, np, pn, pp, nn
## qxy (proportion of x around y): qnn, qpp, qnp, qpn
n, p, np, pn, pp, nn, qnn, qpp, qnp, qpn = var('n,p,np,pn,pp,nn,qnn,qpp,qnp,qpn')
qnp == np / n
 np/n
qnp = np / n
qpn = pn / p
qnn = nn / n
qpp = pp / p

# Clustering
c = var('c')
## Sonia's version:
eq = c == (qnp + qpn + qnn + qpp) / ( n + p)
eq.full_simplify().factor()

## Test:
eq = c == (qnp + qpn + qnn + qpp) / ( (n + p) * (n + p))
eq.full_simplify().factor()

