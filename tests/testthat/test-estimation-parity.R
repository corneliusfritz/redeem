library(testthat)
library(redeem)

test_that("NR and MM estimation methods achieve parity (Simple REM)", {
  set.seed(42)
  n_nodes <- 6
  n_events <- 50

  times <- cumsum(rexp(n_events, rate = 1))
  from <- sample(1:n_nodes, n_events, replace = TRUE)
  to <- sample(1:n_nodes, n_events, replace = TRUE)
  valid <- from != to
  events <- cbind(times[valid], from[valid], to[valid])
  colnames(events) <- c("time", "from", "to")

  fit_nr <- rem(
    formula = ~ 1,
    events = events,
    n_nodes = n_nodes,
    directed = FALSE,
    control = control.redeem(estimate = "NR", verbose = FALSE)
  )

  fit_mm <- rem(
    formula = ~ 1,
    events = events,
    n_nodes = n_nodes,
    directed = FALSE,
    control = control.redeem(estimate = "Blockwise", verbose = FALSE)
  )

  # Access via $model$coefficients
  coef_nr <- fit_nr$model$coefficients
  coef_mm <- fit_mm$model$coefficients
  expect_equal(fit_nr$model$llh, fit_mm$model$llh, tolerance = 1e-4)
  expect_equal(coef_nr, coef_mm, tolerance = 1e-4)
})

test_that("NR and MM parity for directed models (no degrees)", {
  set.seed(42)
  n_nodes <- 6
  n_events <- 50

  times <- cumsum(rexp(n_events, rate = 1))
  from <- sample(1:n_nodes, n_events, replace = TRUE)
  to <- sample(1:n_nodes, n_events, replace = TRUE)
  valid <- from != to
  events <- cbind(times[valid], from[valid], to[valid])
  colnames(events) <- c("time", "from", "to")

  fit_nr <- rem(
    formula = ~ 1,
    events = events,
    n_nodes = n_nodes,
    directed = TRUE,
    control = control.redeem(estimate = "NR", verbose = FALSE)
  )

  fit_mm <- rem(
    formula = ~ 1,
    events = events,
    n_nodes = n_nodes,
    directed = TRUE,
    control = control.redeem(estimate = "Blockwise", verbose = FALSE)
  )
  expect_equal(fit_nr$model$llh, fit_mm$model$llh, tolerance = 1e-4)
  expect_equal(fit_nr$model$coefficients, fit_mm$model$coefficients, tolerance = 1e-4)
})

test_that("NR and MM parity on complex simulated DEM data", {
  set.seed(123)
  n_nodes <- 6

  # Define a model with inertia and degrees
  formula_0_1 <- ~ 1 + inertia() + degrees
  formula_1_0 <- ~ 1

  # Coefficients: Intercept=1.0, inertia=0.5
  coef_0_1 <- c(1.0, 0.5)
  # Degree effects (one per node for undirected)
  coef_degree_0_1 <- rnorm(n_nodes, 0, 0.2)

  # Intercept for termination
  coef_1_0 <- c(-0.5)

  # Simulate 150 events
  sim_events <- dem.simulate(
    formula_0_1 = formula_0_1,
    formula_1_0 = formula_1_0,
    coef_0_1 = coef_0_1,
    coef_1_0 = coef_1_0,
    coef_degree_0_1 = coef_degree_0_1,
    n_events = 150,
    n_nodes = n_nodes,
    directed = FALSE,
    seed = 42,
    verbose = FALSE
  )

  # 1. Newton-Raphson
  fit_nr <- dem(
    formula_0_1 = formula_0_1,
    formula_1_0 = formula_1_0,
    events = sim_events,
    n_nodes = n_nodes,
    directed = FALSE,
    control = control.redeem(
      estimate = "NR",
      it_max = 50,
      tol = 1e-8,
      verbose = FALSE
    )
  )

  # 2. Blockwise MM
  fit_mm <- dem(
    formula_0_1 = formula_0_1,
    formula_1_0 = formula_1_0,
    events = sim_events,
    n_nodes = n_nodes,
    directed = FALSE,
    control = control.redeem(
      estimate = "Blockwise",
      it_max = 200,
      tol = 1e-8,
      accelerated = TRUE,
      verbose = FALSE
    )
  )

  # Access coefficients via the sub-models as requested
  coef_nr <- fit_nr$model_0_1$coefficients
  coef_mm <- fit_mm$model_0_1$coefficients

  common_names <- intersect(names(coef_nr), names(coef_mm))
  expect_gt(length(common_names), 5)

  vals_nr <- coef_nr[common_names]
  vals_mm <- coef_mm[common_names]

  valid <- is.finite(vals_nr) & is.finite(vals_mm)

  # Strict correlation check for parity
  expect_gt(cor(vals_nr[valid], vals_mm[valid]), 0.999)

  # Coefficients should be numerically very close
  expect_equal(vals_nr[valid], vals_mm[valid], tolerance = 1e-4)

  # LLH Parity
  expect_equal(fit_nr$model_0_1$llh, fit_mm$model_0_1$llh, tolerance = 1e-4)
})

test_that("NR and MM parity on complex simulated REM data", {
  set.seed(456)
  n_nodes <- 6

  formula <- ~ 1 + inertia() + degrees

  coef_rem <- c(0.5, 0.3)
  coef_degree <- rnorm(n_nodes, 0, 0.1)

  sim_events <- dem.simulate(
    formula_0_1 = formula,
    coef_0_1 = coef_rem,
    coef_degree_0_1 = coef_degree,
    n_events = 100,
    n_nodes = n_nodes,
    directed = FALSE,
    seed = 123
  )

  fit_nr <- rem(
    formula = formula,
    events = sim_events[, 1:3],
    n_nodes = n_nodes,
    directed = FALSE,
    control = control.redeem(estimate = "NR", verbose = FALSE)
  )

  fit_mm <- rem(
    formula = formula,
    events = sim_events[, 1:3],
    n_nodes = n_nodes,
    directed = FALSE,
    control = control.redeem(estimate = "Blockwise", verbose = FALSE)
  )

  # Access via $model$coefficients as requested
  coef_nr <- fit_nr$model$coefficients
  coef_mm <- fit_mm$model$coefficients

  common_names <- intersect(names(coef_nr), names(coef_mm))
  valid <- is.finite(coef_nr[common_names]) & is.finite(coef_mm[common_names])

  expect_gt(cor(coef_nr[common_names][valid], coef_mm[common_names][valid]), 0.999)
  expect_equal(coef_nr[common_names][valid], coef_mm[common_names][valid], tolerance = 1e-4)

  # LLH Parity
  expect_equal(fit_nr$model$llh, fit_mm$model$llh, tolerance = 1e-4)
})

