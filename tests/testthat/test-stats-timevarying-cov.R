library(testthat)
library(redeem)



test_that("preprocess handles static covariates as before", {
  n <- 5
  events <- matrix(c(
    1.2, 1, 2, 1,
    2.5, 1, 2, 0,
    3.1, 2, 3, 1,
    4.4, 2, 3, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")
  
  cov_mat <- matrix(runif(n*n), n, n)
  
  # Static version
  res_static <- redeem:::preprocess(
    edgelist = events,
    terms = "dyadic_cov",
    data_list = list(cov_mat),
    transformations = "identity",
    n_nodes = n,
    verbose = FALSE,
    directed = FALSE,
    simultaneous_interactions = FALSE
  )
  
  # Time-varying version with single element at t=0
  res_tv <- redeem:::preprocess(
    edgelist = events,
    terms = "dyadic_cov",
    data_list = list(list("0" = cov_mat)),
    transformations = "identity",
    n_nodes = n,
    verbose = FALSE,
    directed = FALSE,
    simultaneous_interactions = FALSE
  )
  
  expect_equal(res_static, res_tv)
})

test_that("preprocess handles time-varying covariates correctly", {
  n <- 3
  events <- matrix(c(
    1, 1, 2, 1,
    5, 1, 2, 0,
    10, 2, 3, 1,
    15, 2, 3, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")
  
  cov0 <- matrix(0, n, n)
  cov10 <- matrix(1, n, n) # Change at t=10
  
  res <- redeem:::preprocess(
    edgelist = events,
    terms = "dyadic_cov",
    data_list = list(list("0" = cov0, "10" = cov10)),
    transformations = "identity",
    n_nodes = n,
    verbose = FALSE,
    directed = FALSE,
    simultaneous_interactions = FALSE
  )
  
  # Check that for events before t=10, the covariate was cov0 (0)
  # Preprocessed matrix columns: time_new, time, pair_id, status, event, from, to, from_avail, to_avail, dyadic_cov
  # dyadic_cov is at index 10 (1-indexed)
  
  # Events are:
  # 1. t=1 (Start). Interval [0, 1]. Stats from 0 (0).
  # 2. t=5 (End). Interval [1, 5]. Stats from 1 (0).
  # 3. t=10 (Start). Interval [5, 10]. Stats from 5 (0).
  # 4. t=15 (End). Interval [10, 15]. Stats from 10 (1).
  
  events_only <- res[res[, 5] == 1, ]
  
  expect_equal(events_only[1, 10], 0) # Event at t=1
  expect_equal(events_only[2, 10], 0) # Event at t=5
  expect_equal(events_only[3, 10], 0) # Event at t=10
  expect_equal(events_only[4, 10], 1) # Event at t=15 (first event seeing cov10)
})

test_that("preprocess throws error if first time is after first event", {
  n <- 3
  events <- matrix(c(
    1, 1, 2, 1,
    5, 1, 2, 0
  ), ncol = 4, byrow = TRUE)
  
  cov10 <- matrix(1, n, n)
  
  expect_error(
    redeem:::preprocess(
      edgelist = events,
      terms = "dyadic_cov",
      data_list = list(list("10" = cov10)),
      transformations = "identity",
      n_nodes = n,
      verbose = FALSE,
      directed = FALSE,
      simultaneous_interactions = FALSE
    ),
    "must be at or before the first event time"
  )
})
