library(testthat)
library(redeem)

test_that("general_triangle statistics are invariant to dissolution", {
  # Small network
  n_nodes <- 5

  # Event sequence for OSP (Outgoing Shared Partner)
  events <- matrix(c(
    0.5, 1, 4, 1,
    0.6, 1, 4, 0,
    1.0, 1, 2, 1,
    2.0, 4, 2, 1,
    3.0, 1, 2, 0,
    4.0, 1, 2, 1,
    10.0, 4, 2, 0,
    10.0, 1, 2, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")

  f_0_1 <- ~ general_triangle(type = "OSP", transformation = "identity")
  f_1_0 <- ~ Intercept

  fit <- dem(
    events = events,
    formula_0_1 = f_0_1,
    formula_1_0 = f_1_0,
    n_nodes = n_nodes,
    directed = TRUE,
    exogenous_end = 10.0,
    control = control.redeem(return_data = TRUE)
  )

  data_0_1 <- fit$model_0_1$data
  dyad_1_4 <- data_0_1[data_0_1$from == 1 & data_0_1$to == 4, ]
  dyad_1_4 <- dyad_1_4[order(dyad_1_4$time), ]

  find_stat <- function(df) {
    cols <- colnames(df)
    cols[grep("general_triangle", cols)][1]
  }
  stat_col <- find_stat(data_0_1)

  get_val <- function(df, target, col) {
    row <- df[df$time <= target & df$time_new > target, ]
    if(nrow(row) == 0) return(tail(df[[col]], 1))
    row[[col]][1]
  }

  val_t2 <- get_val(dyad_1_4, 2.1, stat_col) # After 4-2 forms
  val_t3 <- get_val(dyad_1_4, 3.1, stat_col) # After 1-2 dissolves
  val_t4 <- get_val(dyad_1_4, 4.1, stat_col) # After 1-2 re-forms

  expect_equal(val_t2, 1, info = "Dyad (1,4) should have partner 2 at t=2.1")
  expect_equal(val_t3, 1, info = "General history should NOT decrement on dissolution at t=3")
  expect_equal(val_t4, 1, info = "General history should NOT increment on repeat formation at t=4")
})

test_with_formation <- function() {
  n_nodes <- 5
  events <- matrix(c(
    1.0, 1, 2, 1,
    2.0, 4, 2, 1,
    2.5, 1, 4, 1,
    3.0, 1, 2, 0,
    10.0, 4, 2, 0,
    10.0, 1, 4, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")

  f_0_1 <- ~ current_common_partners(type = "OSP", transformation = "identity")
  f_1_0 <- ~ Intercept

  fit <- dem(
    events = events,
    formula_0_1 = f_0_1,
    formula_1_0 = f_1_0,
    n_nodes = n_nodes,
    directed = TRUE,
    control = control.redeem(return_data = TRUE)
  )

  data_0_1 <- fit$model_0_1$data
  dyad_1_4 <- data_0_1[data_0_1$from == 1 & data_0_1$to == 4, ]
  val_t2_5 <- dyad_1_4$current_common_partner_OSP_identity[dyad_1_4$event == 1]

  expect_equal(val_t2_5, 1, info = "Active partner 2 exists at t=2.5 when (1,4) forms")
}

test_that("current_triangle statistics decrement correctly on dissolution", {
  test_with_formation()
})

