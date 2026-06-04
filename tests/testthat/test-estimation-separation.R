library(testthat)
library(redeem)
library(data.table)

test_that("detect_separation correctly identifies separated covariates in estimate_mm", {
  set.seed(123)
  n_nodes <- 5
  
  # Create a simple event sequence where node 1 NEVER initiates an event
  # and node 1 is the only one with x=1
  events_no_1 <- matrix(c(
    rep(c(2, 3, 4, 1), 5),
    rep(c(3, 4, 5, 1), 5),
    rep(c(4, 5, 2, 1), 5),
    rep(c(5, 2, 3, 1), 5)
  ), ncol = 4, byrow = TRUE)
  colnames(events_no_1) <- c("time", "from", "to", "type")
  events_no_1[,1] <- seq_len(nrow(events_no_1))
  
  # Node 1 has x=1, others have x=0
  x_vec <- c(1, 0, 0, 0, 0)
  
  # Function to pick sender's covariate value
  sender_x <- function(i, j) i
  
  # Estimate model with this covariate
  fit_sep <- rem(
    events = events_no_1,
    formula = ~ monadic_cov(data = x_vec, fun = sender_x),
    n_nodes = n_nodes,
    directed = TRUE,
    control = control.redeem(estimate = "Blockwise", it_max = 5)
  )
  
  # Check if the coefficient for the separated covariate is very small or -Inf
  coefs <- fit_sep$model$est_core
  expect_true(any(coefs < -2) || any(coefs == -Inf))
})

test_that("detect_separation correctly identifies separated covariates in estimate_mmt", {
  set.seed(123)
  n_nodes <- 5
  events_no_1 <- matrix(c(
    rep(c(2, 3, 4, 1), 5),
    rep(c(3, 4, 5, 1), 5),
    rep(c(4, 5, 2, 1), 5),
    rep(c(5, 2, 3, 1), 5)
  ), ncol = 4, byrow = TRUE)
  colnames(events_no_1) <- c("time", "from", "to", "type")
  events_no_1[,1] <- seq_len(nrow(events_no_1))
  
  x_vec <- c(1, 0, 0, 0, 0)
  sender_x <- function(i, j) i
  
  # Estimate model with time-varying baseline to trigger estimate_mmt
  fit_sep_t <- rem(
    events = events_no_1,
    formula = ~ monadic_cov(data = x_vec, fun = sender_x) + baseline(changepoints = 10),
    n_nodes = n_nodes,
    directed = TRUE,
    control = control.redeem(estimate = "Blockwise", it_max = 5)
  )
  
  # Check if the effect is very small or -Inf
  coefs_t <- fit_sep_t$model$est_core
  expect_true(any(coefs_t < -2) || any(coefs_t == -Inf))
})

test_that("inf_unidentifiable = FALSE control option keeps separated/unidentifiable coefficients finite", {
  set.seed(123)
  n_nodes <- 5
  events_no_1 <- matrix(c(
    rep(c(2, 3, 4, 1), 5),
    rep(c(3, 4, 5, 1), 5),
    rep(c(4, 5, 2, 1), 5),
    rep(c(5, 2, 3, 1), 5)
  ), ncol = 4, byrow = TRUE)
  colnames(events_no_1) <- c("time", "from", "to", "type")
  events_no_1[,1] <- seq_len(nrow(events_no_1))

  x_vec <- c(1, 0, 0, 0, 0)
  sender_x <- function(i, j) i

  # Estimate model with inf_unidentifiable = FALSE
  fit_sep_finite <- rem(
    events = events_no_1,
    formula = ~ monadic_cov(data = x_vec, fun = sender_x) + degree + baseline(changepoints = 10),
    n_nodes = n_nodes,
    directed = TRUE,
    control = control.redeem(estimate = "Blockwise", it_max = 5, inf_unidentifiable = FALSE)
  )

  # Check that unidentifiable core coefficients are NOT set to -Inf
  coefs_core <- fit_sep_finite$model$est_core
  expect_false(any(coefs_core == -Inf))

  # Check that zero-count degree estimates are NOT set to -Inf
  coefs_degree <- fit_sep_finite$model$est_degree
  expect_false(any(coefs_degree == -Inf))
})

