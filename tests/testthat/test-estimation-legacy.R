library(testthat)
library(redeem)

test_that("control.redeem supports legacy option", {
  ctrl <- control.redeem(legacy = TRUE)
  expect_true(ctrl$legacy)
  
  ctrl_default <- control.redeem()
  expect_false(ctrl_default$legacy)
})

test_that("dem() works with legacy = TRUE", {
  set.seed(123)
  n_nodes <- 5
  events <- dem.simulate(
    formula_0_1 = ~intercept,
    formula_1_0 = ~intercept,
    coef_0_1 = -2,
    coef_1_0 = -1,
    n_events = 30,
    n_nodes = n_nodes,
    directed = FALSE
  )
  
  # Standard NR
  fit_std <- dem(
    events = events,
    formula_0_1 = ~intercept,
    formula_1_0 = ~intercept,
    n_nodes = n_nodes,
    directed = FALSE,
    control = control.redeem(estimate = "NR", legacy = FALSE)
  )
  
  # Legacy NR
  fit_legacy <- dem(
    events = events,
    formula_0_1 = ~intercept,
    formula_1_0 = ~intercept,
    n_nodes = n_nodes,
    directed = FALSE,
    control = control.redeem(estimate = "NR", legacy = TRUE)
  )
  
  expect_s3_class(fit_legacy, "dem")
  expect_s3_class(fit_legacy$model_0_1, "dem.nr")
  
  # They should be numerically very close for this simple model
  expect_equal(
    fit_std$model_0_1$coefficients,
    fit_legacy$model_0_1$coefficients,
    tolerance = 1e-7
  )
})

test_that("rem() works with legacy = TRUE", {
  set.seed(123)
  n_nodes <- 5
  events <- dem.simulate(
    formula_0_1 = ~intercept,
    formula_1_0 = ~intercept,
    coef_0_1 = -2,
    coef_1_0 = -1,
    n_events = 30,
    n_nodes = n_nodes,
    directed = FALSE
  )
  # Keep only start events for REM
  rem_events <- events[events[, 4] == 1, 1:3]
  
  # Standard NR
  fit_std <- rem(
    events = rem_events,
    formula = ~intercept,
    n_nodes = n_nodes,
    directed = FALSE,
    control = control.redeem(estimate = "NR", legacy = FALSE)
  )
  
  # Legacy NR
  fit_legacy <- rem(
    events = rem_events,
    formula = ~intercept,
    n_nodes = n_nodes,
    directed = FALSE,
    control = control.redeem(estimate = "NR", legacy = TRUE)
  )
  
  expect_s3_class(fit_legacy, "rem")
  expect_s3_class(fit_legacy$model, "dem.nr")
  
  expect_equal(
    fit_std$model$coefficients,
    fit_legacy$model$coefficients,
    tolerance = 1e-7
  )
})
