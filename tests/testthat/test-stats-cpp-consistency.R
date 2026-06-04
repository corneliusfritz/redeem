library(testthat)
library(redeem)
library(data.table)
test_that("Directed Degree terms are consistent between R and C++", {
  n_nodes <- 3
  # Event: 1 -> 2 (forms) at t=1.0. End observation at t=2.0.
  events <- matrix(c(1.0, 1, 2, 1), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")

  fit <- dem(
    events = events,
    formula_0_1 = ~ Intercept + degree_out_sender(),
    n_nodes = n_nodes,
    directed = TRUE,
    exogenous_end = 2.0,
    control = control.redeem(return_data = TRUE)
  )

  dt <- as.data.table(fit$model_0_1$data)

  # After t=1.0, node 1 has out-degree 1.
  # Dyad (1, 3) has node 1 as sender, so it should be included at t=1.0.
  dyad_1_3 <- dt[from == 1 & to == 3 & abs(time - 1.0) < 1e-7]
  expect_true(nrow(dyad_1_3) > 0)
  expect_equal(as.numeric(dyad_1_3$degree_out_sender_identity), 1)
})

test_that("Directed Count terms are consistent", {
  n_nodes <- 3
  events <- matrix(c(
    1.0, 1, 2, 1,
    2.0, 1, 2, 0,
    3.0, 1, 2, 1
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")

  fit <- dem(
    events = events,
    formula_0_1 = ~ Intercept + general_count_out_sender(),
    n_nodes = n_nodes,
    directed = TRUE,
    exogenous_end = 4.0,
    control = control.redeem(return_data = TRUE)
  )

  dt <- as.data.table(fit$model_0_1$data)
  # After t=3.0, node 1 has general out-count 2.
  dyad_1_3_t3 <- dt[from == 1 & to == 3 & abs(time - 3.0) < 1e-7]
  expect_true(nrow(dyad_1_3_t3) > 0)
  expect_equal(as.numeric(dyad_1_3_t3$general_count_out_sender_identity), 2)
})

test_that("Undirected Degree terms are consistent", {
  n_nodes <- 3
  # Event: 1-2 forms at t=1. End at t=2.
  events <- matrix(c(1.0, 1, 2, 1), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")

  fit <- dem(
    events = events,
    formula_0_1 = ~ Intercept + degree_sum(),
    n_nodes = n_nodes,
    directed = FALSE,
    exogenous_end = 2.0,
    control = control.redeem(return_data = TRUE)
  )

  dt <- as.data.table(fit$model_0_1$data)
  # Interval [1.0, 2.0]: deg(1)=1, deg(2)=1, deg(3)=0.
  # Dyad (1, 3) should be included because deg(1) changed.
  expect_equal(as.numeric(dt[from == 1 & to == 3 & abs(time - 1.0) < 1e-7]$degree_sum_identity), 1)
})

test_that("Triadic terms (OTP) are consistent", {
  n_nodes <- 3
  # Sequence: 1->2 (t=1), 2->3 (t=2). End at t=3.
  events <- matrix(c(
    1.0, 1, 2, 1,
    2.0, 2, 3, 1
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")

  fit <- dem(
    events = events,
    formula_0_1 = ~ Intercept + current_common_partners(type = "OTP"),
    n_nodes = n_nodes,
    directed = TRUE,
    exogenous_end = 3.0,
    control = control.redeem(return_data = TRUE)
  )

  dt <- as.data.table(fit$model_0_1$data)
  # Interval [2.0, 3.0]: 1->2 and 2->3 exist. OTP(1,3)=1.
  # Dyad (1, 3) should be included at t=2.0 because 2->3 completes a path for it.
  dyad_1_3 <- dt[from == 1 & to == 3 & abs(time - 2.0) < 1e-7]
  expect_true(nrow(dyad_1_3) > 0)
  expect_equal(as.numeric(dyad_1_3$current_common_partner_OTP_identity), 1)
})

test_that("Structural terms (Inertia, Reciprocity) are consistent", {
  n_nodes <- 2
  # 1->2 forms (t=1), dissolves (t=2), forms (t=3). End at t=4.
  events <- matrix(c(
    1.0, 1, 2, 1,
    2.0, 1, 2, 0,
    3.0, 1, 2, 1
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")

  fit <- dem(
    events = events,
    formula_0_1 = ~ Intercept + inertia() + reciprocity(),
    n_nodes = n_nodes,
    directed = TRUE,
    exogenous_end = 4.0,
    control = control.redeem(return_data = TRUE)
  )

  dt <- as.data.table(fit$model_0_1$data)
  # Dyad (1, 2) is always included when it is the event pair.
  # At t=2.0 (after dissolution), inertia for (1,2) should be 1.
  row2_12 <- dt[from == 1 & to == 2 & abs(time - 2.0) < 1e-7]
  expect_equal(as.numeric(row2_12$inertia_identity), 1)

  # reciprocity(2,1) should be 1 after 1->2 dissolves (t=2.0).
  # In DEM, non-windowed structural terms wait for dissolution.
  row2_21 <- dt[from == 2 & to == 1 & abs(time - 2.0) < 1e-7]
  expect_true(nrow(row2_21) > 0)
  expect_equal(as.numeric(row2_21$reciprocity_identity), 1)
})

test_that("Covariate terms (Dyadic, Monadic) are consistent", {
  n_nodes <- 3
  # Add a type 3 event at t=0 to force all rows at initialization
  events <- matrix(c(
    0.0, 1, 2, 3,
    1.0, 1, 2, 1
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")

  W <- matrix(0, 3, 3)
  W[1, 3] <- 5
  X <- c(10, 20, 30)

  fit <- dem(
    events = events,
    formula_0_1 = ~ Intercept + dyadic_cov(W) + monadic_cov(X, fun = function(s, r) s),
    n_nodes = n_nodes,
    directed = TRUE,
    exogenous_end = 2.0,
    control = control.redeem(return_data = TRUE)
  )

  dt <- as.data.table(fit$model_0_1$data)

  # Covariates are verified by absence of errors in fit and manual check of column names.
  # Static covariates are consistent if the model matrix can be formed correctly.
  expect_true("dyadic_cov_W_identity" %in% names(dt))
  expect_true("monadic_cov_X_identity" %in% names(dt))
})
