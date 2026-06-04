library(testthat)
library(redeem)

test_that("Constant baseline simulation works", {
  n_nodes <- 5
  # 1. Using Intercept in coef
  set.seed(123)
  sim1 <- dem.simulate(
    n_nodes = n_nodes, n_events = 50,
    formula_0_1 = ~1, coef_0_1 = c(Intercept = -2)
  )

  # 2. Using baseline argument (now it is empty for constant because Intercept is present)
  set.seed(123)
  sim2 <- dem.simulate(
    n_nodes = n_nodes, n_events = 50,
    formula_0_1 = ~1, baseline_0_1 = numeric(0)
  )
  # Update: if we want to change the intercept via baseline_0_1,
  # and it's full_baseline=FALSE, we must provide changepoints.
  # For constant simulation with non-zero intercept, one should use coef_0_1.

  expect_equal(nrow(sim1), 50)
  expect_equal(nrow(sim2), 50)
})

test_that("Time-varying baseline simulation works", {
  n_nodes <- 5
  tc <- c(5, 10)
  bl <- c(-2, -5) # After 5 is -2, after 10 is -5. Slice [0, 5) is 0.
  # Formula ~1 has intercept, so full_baseline = FALSE.
  # baseline_0_1 length must match tc.

  set.seed(123)
  sim <- dem.simulate(
    n_nodes = n_nodes, time = 15,
    formula_0_1 = ~1 + baseline(changepoints = tc),
    baseline_0_1 = bl
  )

  # Verify that events in interval [10, 15] have lower intensity (less events)
  # LPs:
  # [0, 5]: Intercept=0 + baseline=0 = 0
  # [5, 10]: Intercept=0 + baseline=-2 = -2
  # [10, 15]: Intercept=0 + baseline=-5 = -5

  events_start <- sim[sim[, 1] <= 5, ]
  events_mid <- sim[sim[, 1] > 5 & sim[, 1] <= 10, ]
  events_end <- sim[sim[, 1] > 10, ]

  rate_start <- nrow(events_start) / 5
  rate_mid <- nrow(events_mid) / 5
  rate_end <- nrow(events_end) / 5

  expect_true(rate_start > rate_mid)
  expect_true(rate_mid > rate_end)
})

test_that("Constant baseline in REM simulation works", {
  n_nodes <- 5
  set.seed(123)
  sim1 <- rem.simulate(
    n_nodes = n_nodes, n_events = 50,
    formula = ~1, coef = c(Intercept = -2)
  )

  set.seed(123)
  sim2 <- rem.simulate(
    n_nodes = n_nodes, n_events = 50,
    formula = ~1, baseline = numeric(0)
  )

  expect_equal(nrow(sim1), 50)
  expect_equal(nrow(sim2), 50)
})

test_that("Complex baseline + covariate changepoints work", {
  n_nodes <- 5

  # Covariate changes at 4 and 8
  cov_data <- list(
    "0" = matrix(1, n_nodes, n_nodes),
    "4" = matrix(2, n_nodes, n_nodes),
    "8" = matrix(1, n_nodes, n_nodes)
  )

  # Baseline changes at 6
  tc <- 6
  bl <- c(-5) # After 6

  # Formula ~ cov has intercept by default in dem.simulate if simulation=TRUE
  set.seed(123)
  sim <- dem.simulate(
    n_nodes = n_nodes, time = 12,
    formula_0_1 = ~ intercept + dyadic_cov(cov_data) + baseline(changepoints = tc),
    baseline_0_1 = bl
  )

  # LPs:
  # [0, 4]: Intercept=0, cov=1, bl=0 -> LP = 0 + 1*1.5 + 0 = 1.5 (if coef=1.5)
  # [4, 6]: Intercept=0, cov=2, bl=0 -> LP = 0 + 2*1.5 + 0 = 3.0
  # [6, 8]: Intercept=0, cov=2, bl=-5 -> LP = 0 + 2*1.5 - 5 = -2.0
  # [8, 12]: Intercept=0, cov=1, bl=-5 -> LP = 0 + 1*1.5 - 5 = -3.5

  set.seed(123)
  sim_c <- dem.simulate(
    n_nodes = n_nodes, time = 12,
    formula_0_1 = ~ intercept + dyadic_cov(cov_data) + baseline(changepoints = tc),
    coef_0_1 = c("dyadic_cov(cov_data)" = 1.5),
    baseline_0_1 = bl
  )

  rate_1 <- nrow(sim_c[sim_c[, 1] <= 4, ]) / 4
  rate_2 <- nrow(sim_c[sim_c[, 1] > 4 & sim_c[, 1] <= 6, ]) / 2
  rate_3 <- nrow(sim_c[sim_c[, 1] > 6 & sim_c[, 1] <= 8, ]) / 2
  rate_4 <- nrow(sim_c[sim_c[, 1] > 8, ]) / 4

  expect_true(rate_2 > rate_1)
  expect_true(rate_2 > rate_3)
  expect_true(rate_3 > rate_4)
})
