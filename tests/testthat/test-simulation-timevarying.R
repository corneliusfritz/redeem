library(testthat)
library(redeem)

# Time-varying Simulation and Prediction tests

test_that("rem.simulate supports time-varying covariates", {
  n_nodes <- 20
  # 1. Setup time-varying covariates
  set.seed(42)
  # X1: low values, X2: high values
  X1 <- matrix(runif(n_nodes^2, 0, 1), n_nodes, n_nodes)
  X2 <- matrix(runif(n_nodes^2, 5, 10), n_nodes, n_nodes)
  X_tv <- list("0" = X1, "0.5" = X2)
  
  formula_tv <- ~ intercept + dyadic_cov(data = X_tv)
  coef <- c(-2, 1) # Intercept -2, Covariate 1. Total rate will jump at 0.5.
  
  events <- rem.simulate(formula = formula_tv, coef = coef, n_nodes = n_nodes, time = 1.0, seed = 123)
  
  expect_true(nrow(events) > 0)
  
  # Verify that we have significantly more events in the second half [0.5, 1.0]
  events_post <- events[events[, 1] >= 0.5, ]
  
  # We use a proportion test: we expect > 80% of events in the second half
  expect_gt(nrow(events_post) / nrow(events), 0.8)
})

test_that("dem.simulate supports time-varying covariates", {
  n_nodes <- 10
  X_low <- matrix(0, n_nodes, n_nodes)
  X_high <- matrix(10, n_nodes, n_nodes)
  X_tv <- list("0.2" = X_high, "0.4" = X_low, "0" = X_low)
  
  formula_tv <- ~ intercept + dyadic_cov(data = X_tv)
  
  # Simulate DEM
  # Intercept low, covariate high
  events <- dem.simulate(formula_0_1 = formula_tv, formula_1_0 = ~ intercept,
                        coef_0_1 = c(-5, 1), coef_1_0 = c(0),
                        n_nodes = n_nodes, time = 0.5, seed = 1)
  
  expect_true(nrow(events) > 0)
  
  # Check if events mostly happen in [0.2, 0.4] where X_high is active
  events_active <- events[events[, 1] >= 0.2 & events[, 1] < 0.4, ]
  
  # We expect the vast majority of events to happen in the active window
  # Given the massive rate difference (exp(5) vs exp(-5)), the ratio should be safely above 0.8,
  # but some termination events might naturally fall outside the window.
  expect_gt(nrow(events_active) / nrow(events), 0.8)
})

test_that("get_ranking handles time-varying covariates (prediction/validation)", {
  n_nodes <- 5
  X_0 <- matrix(0, n_nodes, n_nodes)
  X_1 <- matrix(1, n_nodes, n_nodes)
  X_tv <- list("0" = X_0, "1.0" = X_1)
  
  formula_rem <- ~ intercept + dyadic_cov(data = X_tv)
  coef <- c(-2, 1)
  
  edgelist_train <- matrix(c(0.1, 1, 2, 0,
                            0.2, 2, 3, 0), ncol = 4, byrow = TRUE)
  
  edgelist_test <- matrix(c(0.6, 1, 3,
                           0.7, 4, 5), ncol = 3, byrow = TRUE)
  
  mock_rem <- list(
    formula = formula_rem,
    model = list(coef = coef),
    n_nodes = n_nodes,
    events = edgelist_train,
    directed = FALSE
  )
  class(mock_rem) <- "rem"
  
  res <- get_ranking(
    object = mock_rem,
    edgelist_test = edgelist_test,
    k_max = 10
  )
  
  expect_s3_class(res, "ranking_redeem")
})

test_that("safe_stod throws descriptive error for non-numeric keys", {
  n_nodes <- 5
  X_tv <- list("zero" = matrix(0, n_nodes, n_nodes))
  
  expect_error(
    rem.simulate(formula = ~ dyadic_cov(data = X_tv), coef = 1, n_nodes = n_nodes, time = 1),
    "All names of the 'data' list must be numeric"
  )
})
