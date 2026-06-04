library(testthat)
library(redeem)

test_that("current_interaction advances at covariate-only changepoints", {
  n_nodes <- 3
  # Event at t=0: interaction starts
  # Covariate change at t=5: no event
  # Event at t=10: interaction ends
  
  tv_data <- list("0" = matrix(0, 3, 3), "5" = matrix(1, 3, 3))
  
  # edgelist: time, from, to, type
  edgelist <- matrix(c(
    0, 1, 2, 1,   # interaction starts (type 1)
    10, 1, 2, 0   # interaction ends (type 0)
  ), 2, 4, byrow=TRUE)
  
  # Pass time-varying intercept to force changepoints into the engine
  # We make it CHANGE at t=5 so it's not merged!
  intercept_tv <- list("0" = matrix(1, 3, 3), "5" = matrix(1.1, 3, 3))
  
  res <- redeem:::preprocess(edgelist, 
                    terms = c("intercept", "current_interaction"), 
                    data_list = list(intercept_tv, matrix(0, 3, 3)), 
                    transformations = c("identity", "identity"),
                    n_nodes = n_nodes, verbose = FALSE, directed = TRUE,
                    simultaneous_interactions = TRUE)
  
  # Column names: time_new, time, i, j, pair_id, from, to, from_avail, to_avail, intercept, current_interaction
  # current_interaction should be at index 11 (1-indexed)
  
  # Snapshot at t=0 (interaction starts)
  rows_at_0 <- res[res[, 2] == 0, ]
  pair12_at_0 <- rows_at_0[rows_at_0[, 6] == 1 & rows_at_0[, 7] == 2, , drop = FALSE]
  expect_equal(pair12_at_0[1, 11], 0)
  
  # Snapshot at t=5 (covariate changepoint)
  # The duration was advanced by 5 at the beginning of t=5 iteration.
  rows_at_5 <- res[res[, 2] == 5, ]
  pair12_at_5 <- rows_at_5[rows_at_5[, 6] == 1 & rows_at_5[, 7] == 2, , drop = FALSE]
  
  # Duration at t=5 (reflecting interval [5, 10]) should be 5
  expect_equal(pair12_at_5[1, 11], 5)
})

test_that("all intermediate changepoints trigger snapshots", {
  n_nodes <- 3
  # One event at 0, one at 10. Covariate at 3, 7.
  # Intermediate snapshots at 3 and 7 should exist.
  # We MUST make the covariate change at these times, otherwise they are merged.
  intercept_tv <- list("0"=matrix(1,3,3), "3"=matrix(1.1,3,3), "7"=matrix(1.2,3,3))
  edgelist <- matrix(c(0,1,2,1, 10,2,3,1), 2, 4, byrow=TRUE)
  
  res <- redeem:::preprocess(edgelist, terms = "intercept", data_list = list(intercept_tv), 
                    transformations = "identity", n_nodes = n_nodes, verbose = FALSE, directed = TRUE,
                    simultaneous_interactions = TRUE)
  
  unique_times <- unique(res[, 2])
  # We expect 0, 3, 7. 10 is the very last point and has no interval starting there.
  expect_true(all(c(0, 3, 7) %in% unique_times))
})
