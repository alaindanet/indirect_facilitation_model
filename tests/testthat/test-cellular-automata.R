context("cellular-automata")

test_that("cellular automata runs", {
  mod <- ca_two_facilitation_model()
  times(mod)["to"] <- c(to = 3)
  mod2 <- mod 
  solver(mod) <- "myiteration"
  solver(mod2) <- "myiteration2"
  parms(mod)
  mod_run <- sim(mod)
  expect_is(mod_run, "gridModel")
  # Correct time
  expect_equal(out(mod_run) %>% nrow(.), 4)

  mod_run <- sim(mod2)
  expect_is(mod_run, "gridModel")
  # Correct time
  expect_equal(out(mod_run) %>% nrow(.), 4)
  
})
