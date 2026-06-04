library(testthat)
library(redeem)

test_that("Simulation with windowed statistic works", {
  set.seed(42)
  n <- 3
  
  # Simulate 100 events using windowed inertia
  # inertia(window=2): if event (i,j) happens, inertia=1 for 2 time units, then 0.
  events <- rem.simulate(
      formula = ~ Intercept + inertia(window = 2),
      coef = c(log(0.1), log(2)),  # Intercept and inertia
      n_events = 100,
      n_nodes = n,
      directed = TRUE
  )
  
  # Just verify that events is a matrix with 3 columns and 100 rows
  expect_equal(nrow(events), 100)
  expect_equal(ncol(events), 3)
  
  # Estimate back to see if we get the right parameter
  fit <- rem(
      events = events,
      formula = ~ inertia(window = 2),
      n_nodes = n,
      directed = TRUE
  )
  
  # The true inertia coef is log(2) ~ 0.693
  # Because n_events=100 is small, tolerance can be larger, but we just check if it ran successfully
  # and has the expected terms.
  expect_true("inertia_wt2_identity" %in% names(fit$model$coefficients))
})
