context("running_methods")

library(tibble)
library(magrittr)
set.seed(123)

run <- tibble(
  time = seq(0, 200),
  N = runif(n = 201),
  P = runif(n = 201)
  )
bad_run <- run %>%
  dplyr::mutate(
    N = c(NaN, runif(n = 200)),
    P = c(-400, runif(n = 200))
    )
test_that("is_normal_run works well", {
  expect_false(is_run_normal(bad_run))
  expect_true(is_run_normal(run))
    })
test_that("avg_runs return a good summary of the simulation", {

my_calc <- dplyr::slice(run, (n() - 100) :n()) %>%
  tidyr::gather(species, rho, -time) %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(avg = mean(rho)) %>%
  tidyr::spread(species, avg)

  summary_run <- avg_runs(run, cut_row = 100)
  expect_equal(summary_run$N, my_calc$N, tolerance = .001)
  expect_equal(summary_run$P, my_calc$P, tolerance = .001)
  expect_true(summary_run$status)
  # Bad run:
  summary_run <- avg_runs(bad_run)
  expect_false(summary_run$status)
    })
test_that("avg_run returns NA if simulations has stoped early", {
  short_run <- tibble(
    time = seq(0, 2),
    N = runif(n = 3),
    P = runif(n = 3)
    )
  expect_true(avg_runs(short_run)$status)

  short_bad_run <- tibble(
    time = seq(0, 2),
    N = rep(NaN, n = 3),
    P = rep(NaN, n = 3)
    )
  expect_false(avg_runs(short_bad_run)$status)
  })
test_that("the good states are returned", {
  expect_match(def_state(nurse = .3, protegee = .3, sim_status = TRUE),
    "coexistence")
  expect_match(def_state(nurse = .3, protegee = .3, sim_status = FALSE),
    "warning")
  expect_match(def_state(nurse = 3.84e-06, protegee = .3, sim_status = TRUE),
    "protegee")
  expect_match(def_state(nurse = .001, protegee = .001, sim_status = TRUE),
    "desert")
  expect_match(def_state(nurse = .01, protegee = .001, sim_status = TRUE),
    "nurse")
  })

two_d_sim <- run_2d_gradient(gradienty = 0.1,
  time_seq = c(from = 0, to = 1, by = 1))

test_that("Run_2d_gradient returns a list", {
  expect_is(two_d_sim, "list")
  expect_is(two_d_sim[["param"]], "list")
  expect_output(str(two_d_sim), "List of 3")
  })

avg_two_d_sim <- avg_runs(two_d_sim, cut_row = 1)
test_that("Avg_runs returns a list", {

  expect_is(avg_two_d_sim, "list")
  expect_is(avg_two_d_sim[["param"]], "list")
  })
# Bifurcation
test_that("run_bifurc_model runs", {
bifurc_sim <- run_bifurc_model(
  x = 0.1, scenario = "together",
  name_x = "b",
  model = indirect_facilitation_model()
  )
expect_output(str(bifurc_sim), "data.frame")

  })
#Scenarii
test_that("scenarii are well specified", {
  expect_output(str(init_scenarii()), "List of 1")
  expect_output(str(init_scenarii(type = "all")), "List of 6")
  expect_is(names(init_scenarii(type = "all")), "character")

  })
test_that("run_simecol works", {
  gradient <- list(b = c(.8,.2), c = .2)
  scenarii <- init_scenarii(type = "together")
  gradient$scenario <- names(scenarii)
  scenar_gradient <- gradient

  comb <- expand.grid(scenar_gradient) %>%
    dplyr::mutate(
      inits = purrr::map(scenario, function(x) scenarii[[x]])
      )
  param_combination <- dplyr::select(comb, -scenario, -inits) %>% df2list(.)

  run <- Map(run_simecol,
    inits = comb[["inits"]],
    param = param_combination,
    MoreArgs = list(model = two_facilitation_model())
    )
  # Problem in the full function, at Map(), see print()
  expect_is(run, "list")
  expect_length(run, 2)
  })

  u0 <- run_scenarii_gradient(
    gradient = list(g = c(0, 0.1), u = 0),
    model_spec = "two_facilitation_model",
    time_seq = c(from = 0, to = 1, by = 1),
    solver_type = NULL
    )
test_that("run_scenarii_gradient run", {
  expect_is(u0, "scenarii")
  expect_length(u0, 5)

  u_tail <- run_scenarii_gradient(
    gradient = list(g = c(0), b = c(.5, .8)),
    model_spec = "two_facilitation_model",
    time_seq = c(from = 0, to = 10, by = 1),
    solver_type = NULL, set_tail = 5
    )
  expect_equal(nrow(u_tail$run$run[[1]]), 5)

  u_tail_avg <- avg_runs(u_tail, cut_row = 5)
  expect_equal(nrow(u_tail_avg$run), 2)
  expect_equal(which(is.na(u_tail_avg$run)) %>% length(.), 0)

  })
test_that("the filter of runs works", {
  u0_filtered <- filter(u0, g == 0)

  expect_is(u0_filtered, "scenarii")
  expect_equal(unique(u0_filtered$run$g), 0)

  expect_is(u0_filtered$gradient, "list")
  expect_equal(u0_filtered$gradient$g, 0)
  expect_length(u0_filtered, 5)

  u0_avg <- avg_runs(u0)
  u0_avg_filtered <- filter(u0_avg, g == 0)
  expect_is(u0_filtered$gradient, "list")
  expect_equal(u0_filtered$gradient$g, 0)
  expect_length(u0_filtered, 5)

  })

u5 <- run_scenarii_gradient(
  gradient = list(g = c(0, 0.1), u = 5),
  model_spec = "two_facilitation_model",
  time_seq = c(from = 0, to = 1, by = 1),
  solver_type = NULL
  )

test_that("binding scenarii works", {
  u_miss <- u0; u_miss$param["c"]  <- .3
  expect_warning(bind_scenar(u0, u_miss))
  binded <- bind_scenar(u5, u0)
  expect_output(class(binded), "scenarii")
  })
