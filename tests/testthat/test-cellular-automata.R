context("cellular-automata")

test_that("cellular automata runs", {
  mod <- ca_two_facilitation_model()
  times(mod)["to"] <- c(to = 3)
  solver(mod) <- "ca_solver"
  mod_run <- sim(mod)
  expect_is(mod_run, "gridModel")
  # Correct time
  expect_equal(out(mod_run) %>% nrow(.), 4)

})


test_that("run_scenarii_gradient runs", {
  u0 <- run_scenarii_gradient(
    gradient = list(g = c(0, 0.1), b = c(.5, .8)),
    model_spec = "ca_two_facilitation_model",
    time_seq = c(from = 0, to = 1, by = 1),
    solver_type = NULL
    )
})
