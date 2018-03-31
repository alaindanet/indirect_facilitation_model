
devtools::document()
library(simecol)

upca <- upca_model()
parms(upca)
equations(upca)$f <- equations(upca)$f2
test <- sim(upca)
plotupca(test)

# Test
devtools::use_testthat()
testthat

library(devtools)
version
