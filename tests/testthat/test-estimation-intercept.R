library(testthat)
library(redeem)

test_that("Simulation: Intercept is NOT added when baseline() is missing", {
  set.seed(123)
  n_nodes <- 5
  
  # formula without baseline
  formula <- ~ inertia
  
  # We use the internal formula_preprocess_single to check the results
  pre <- redeem:::formula_preprocess_single(formula, n_nodes = n_nodes, simulation = TRUE)
  
  # Check if Intercept is in coef_names
  expect_false("Intercept" %in% names(pre$coef_names))
  expect_true("inertia" %in% names(pre$coef_names))
})

test_that("Simulation: Intercept IS added when baseline() is present", {
  set.seed(123)
  n_nodes <- 5
  
  # formula with baseline
  formula <- ~ baseline(changepoints = c(10)) + inertia
  
  pre <- redeem:::formula_preprocess_single(formula, n_nodes = n_nodes, simulation = TRUE)
  
  # Check if Intercept is in coef_names
  expect_true("Intercept" %in% names(pre$coef_names))
  expect_true("inertia" %in% names(pre$coef_names))
})

test_that("Estimation: Intercept is NOT added by default even if degrees are missing", {
  set.seed(123)
  n_nodes <- 5
  
  # formula without degrees and without baseline
  formula <- ~ inertia
  
  pre <- redeem:::formula_preprocess_single(formula, n_nodes = n_nodes, simulation = FALSE)
  
  # Check if Intercept is in coef_names
  expect_false("Intercept" %in% names(pre$coef_names))
})

test_that("Estimation: Intercept is NOT added when baseline() is present", {
  set.seed(123)
  n_nodes <- 5
  
  # formula with baseline
  formula <- ~ baseline(changepoints = c(10)) + inertia
  
  pre <- redeem:::formula_preprocess_single(formula, n_nodes = n_nodes, simulation = FALSE)
  
  # Check if Intercept is NOT in coef_names (per latest user requirement)
  expect_false("Intercept" %in% names(pre$coef_names))
})

test_that("Simulation: Intercept IS added when TV covariates are present", {
  set.seed(123)
  n_nodes <- 5
  
  # formula with TV covariate (must include time 0)
  tv_data <- list("0" = matrix(0, n_nodes, n_nodes), "10" = matrix(1, n_nodes, n_nodes))
  formula <- ~ dyadic_cov(data = tv_data)
  
  pre <- redeem:::formula_preprocess_single(formula, n_nodes = n_nodes, simulation = TRUE)
  
  # Check if Intercept is in coef_names
  expect_true("Intercept" %in% names(pre$coef_names))
})
