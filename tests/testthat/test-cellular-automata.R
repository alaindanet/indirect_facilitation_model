context("cellular-automata")

test_that("cellular automata runs", {
  mod <- ca_two_facilitation_model()
  solver(mod) <- "myiteration"
  times(mod)["to"] <- c(to = 3)
  parms(mod)
  mod_run <- sim(mod)
  expect_is(mod_run, "gridModel")
  
})
