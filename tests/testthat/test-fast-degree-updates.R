library(testthat)
library(redeem)

test_that("eval_llh_pois handles zero means and outcomes correctly", {
  outcome <- c(1, 0, 1, 0)
  mean <- c(1, 1, 0, 0)
  
  # Term 1: 1*log(1) - 1 = -1
  # Term 2: 0*log(1) - 1 = -1
  # Term 3: 1*log(0) - 0 = -Inf
  # Term 4: 0*log(0) - 0 = 0 (handled as 0 in our loop)
  # Total: -1 - 1 - Inf + 0 = -Inf
  expect_equal(redeem:::eval_llh_pois(outcome, mean, numeric(0)), -Inf)
  
  # Case where it's finite:
  outcome_f <- c(1, 0)
  mean_f <- c(2, 2)
  # 1*log(2) - 2 + 0*log(2) - 2 = log(2) - 4
  expect_equal(redeem:::eval_llh_pois(outcome_f, mean_f, numeric(0)), log(2) - 4)
  
  # Weighted case
  weights <- c(2, 3)
  # 2*(1*log(2) - 2) + 3*(0*log(2) - 2) = 2*log(2) - 4 - 6 = 2*log(2) - 10
  expect_equal(redeem:::eval_llh_pois(outcome_f, mean_f, weights), 2*log(2) - 10)
})

test_that("eval_llh_pois_log handles -Inf log_means correctly", {
  outcome <- c(1, 0, 1, 0)
  log_mean <- c(0, 0, -Inf, -Inf)
  
  # Term 1: 1*0 - exp(0) = -1
  # Term 2: 0*0 - exp(0) = -1
  # Term 3: 1*(-Inf) - exp(-Inf) = -Inf
  # Term 4: 0*(-Inf) - exp(-Inf) = 0
  # Total: -Inf
  expect_equal(redeem:::eval_llh_pois_log(outcome, log_mean, numeric(0)), -Inf)
  
  # Finite case
  outcome_f <- c(1, 0)
  log_mean_f <- c(log(2), log(2))
  expect_equal(redeem:::eval_llh_pois_log(outcome_f, log_mean_f, numeric(0)), log(2) - 4)
})

test_that("update_degree_fast parity and edge cases", {
  set.seed(123)
  n_nodes <- 4
  from <- c(1, 2, 3, 1)
  to <- c(2, 3, 4, 3)
  event <- c(1, 0, 1, 0)
  prediction <- c(0.5, 0.5, 0.5, 0.5)
  weights <- rep(1, length(from))
  
  # Directed Sender Update
  est_degree_dir <- rep(0, 2 * n_nodes)
  res_fast <- redeem:::update_degree_fast(from, to, event, prediction, weights, est_degree_dir, n_nodes, directed = TRUE, update_sender = TRUE)$est_degree
  
  # Manual calculation for node 1:
  # obs_sum = event[1] + event[4] = 1 + 0 = 1
  # pred_sum = prediction[1] + prediction[4] = 0.5 + 0.5 = 1.0
  # est_degree[1] = 0 + log(1/1) = 0
  expect_equal(res_fast[1], 0)
  
  # Manual calculation for node 2:
  # obs_sum = event[2] = 0
  # pred_sum = prediction[2] = 0.5
  # est_degree[2] = 0 + log(1e-15 / 0.5) = log(2e-15)
  expect_equal(res_fast[2], log(2e-15), tolerance = 1e-10)
  
  # Undirected Update
  est_degree_undir <- rep(0, n_nodes)
  res_undir <- redeem:::update_degree_fast(from, to, event, prediction, weights, est_degree_undir, n_nodes, directed = FALSE)$est_degree
  
  # For node 1 (undirected):
  # obs_sum = event[1] (1,2) + event[4] (1,3) = 1 + 0 = 1
  # pred_sum = prediction[1] + prediction[4] = 0.5 + 0.5 = 1.0
  # est_degree[1] = log(sqrt(exp(0) * 1 / (1/exp(0)))) = log(1) = 0
  expect_equal(res_undir[1], 0)
  
  # Node 4 (undirected):
  # obs_sum = event[3] (3,4) = 1
  # pred_sum = prediction[3] = 0.5
  # est_degree[4] = log(sqrt(exp(0) * 1 / (0.5/exp(0)))) = log(sqrt(2)) = 0.5*log(2)
  expect_equal(res_undir[4], 0.5*log(2))

  # Bounds check
  expect_error(redeem:::update_degree_fast(c(1, 5), c(2, 3), c(1, 1), c(1, 1), numeric(0), rep(0, 4), 4, FALSE))
})

test_that("get_A_B_C_D_E_F_exact bounds checking", {
  n_nodes <- 5
  from <- c(1, 6) # 6 is out of bounds
  to <- c(2, 3)
  weight <- c(1, 1)
  cov <- matrix(0, 2, 0)
  time <- c(1, 1)
  
  expect_error(redeem:::get_A_B_C_D_E_F_exact(from, to, weight, cov, time, 1, n_nodes, TRUE))
  
  # Length mismatch
  expect_error(redeem:::get_A_B_C_D_E_F_exact(c(1,2), c(2), weight, cov, time, 1, n_nodes, TRUE))
})
