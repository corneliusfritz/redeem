library(testthat)
library(redeem)

test_that("preprocess handles empty edgelist gracefully", {
  n_nodes <- 5
  # Empty edgelist: 0x4 matrix
  edgelist <- matrix(0, 0, 4)

  # Static case
  res <- redeem:::preprocess(edgelist, terms = "intercept", data_list = list(matrix(0, 5, 5)),
                    transformations = "identity", n_nodes = n_nodes, verbose = FALSE, directed = TRUE,
                    simultaneous_interactions = TRUE)
  expect_equal(nrow(res), 0)

  # Time-varying case
  tv_data <- list("0" = matrix(0, 5, 5), "10" = matrix(1, 5, 5))
  res_tv <- redeem:::preprocess(edgelist, terms = "dyadic_cov", data_list = list(tv_data),
                       transformations = "identity", n_nodes = n_nodes, verbose = FALSE, directed = TRUE,
                       simultaneous_interactions = TRUE)
  expect_equal(nrow(res_tv), 0)
})

test_that("dyadic_cov preserves values through time-varying snapshots", {
  n_nodes <- 3
  # tv_data defines the covariate values at Different times
  K <- 2
  tv_data <- list(
    "0" = matrix(0, n_nodes, n_nodes),
    "0.5" = matrix(10 / (10 + 2), n_nodes, n_nodes)
  )

  # We use a dummy event at t=1 so we get a snapshot for the interval ending at 1
  # The interval [0.5, 1] should use the values from t=0.5

  edgelist <- matrix(c(1.0, 1, 2, 1), 1, 4)

  res <- redeem:::preprocess(edgelist, terms = "dyadic_cov", data_list = list(tv_data),
                    transformations = "identity", n_nodes = n_nodes, verbose = FALSE, directed = TRUE,
                    simultaneous_interactions = TRUE)

  # Row 1: t=0 (Initial)
  # Row 2: t=0.5 (Covariate change)
  # Row 3: t=1.0 (Interaction start)

  # Interval [0.5, 1.0]: stats should be 0.8333
  rows_at_1 <- res[res[, 2] == 1.0, ]
  expect_true(all(abs(rows_at_1[, 10] - (10/12)) < 1e-6))
})

test_that("rem.simulate advances time correctly through empty covariate intervals", {
  n_nodes <- 3
  # Covariate changes at 0, 10, 20
  tv_data <- list(
    "0" = matrix(-10, n_nodes, n_nodes), # Very negative
    "10" = matrix(-10, n_nodes, n_nodes),
    "20" = matrix(10, n_nodes, n_nodes)  # Positive
  )

  # High intensity at the very end to ensure we get events there
  formula_rem <- ~ intercept + dyadic_cov(data = tv_data)
  coef <- c(-5, 1) # Total intensity: -15 (at t<20), +5 (at t>=20)

  # Simulate up to time 25.
  set.seed(42)
  events <- rem.simulate(formula = formula_rem, coef = coef, n_nodes = n_nodes, time = 25, verbose = FALSE)

  if (nrow(events) > 0) {
    # If we got events, they should mostly be at t > 20
    expect_true(any(events[, 1] > 20))
  } else {
    # If we got no events even at t > 20 with intensity > 0, that's a problem for the test seed
    # But usually intensity 5 (exp(5) events per unit time) should produce many events
    fail("No events simulated despite high intensity at [20, 25]")
  }
})

test_that("sentinels handle negative covariate times correctly", {
  n_nodes <- 3
  # Covariate measurements at -10 and -5
  tv_data <- list(
    "-10" = matrix(1, n_nodes, n_nodes),
    "-5"  = matrix(2, n_nodes, n_nodes)
  )

  # At t=0, it should pick the one at -5 (max t <= 0)
  edgelist <- matrix(c(1.0, 1, 2, 1), 1, 4)
  res <- redeem:::preprocess(edgelist, terms = "dyadic_cov", data_list = list(tv_data),
                    transformations = "identity", n_nodes = n_nodes, verbose = FALSE, directed = TRUE,
                    simultaneous_interactions = TRUE)

  # Column 2 has the time. We check a row where time is -5.
  # The value should be 2.
  val_at_neg5 <- res[res[, 2] == -5, 10][1]
  expect_equal(val_at_neg5, 2)
})

test_that("rem.simulate respects n_events even in time-varying intervals", {
  n_nodes <- 5
  tv_data <- list("0" = matrix(0, 5, 5), "100" = matrix(1, 5, 5))
  formula_rem <- ~ intercept + dyadic_cov(data = tv_data)

  # Request exactly 5 events. The interval [0, 100] is long enough to produce many events.
  set.seed(42)
  events <- rem.simulate(formula = formula_rem, coef = c(0, 0), n_nodes = n_nodes,
                         n_events = 5, verbose = FALSE)

  expect_equal(nrow(events), 5)
})
