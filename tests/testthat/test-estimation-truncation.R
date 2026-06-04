
library(testthat)
library(redeem)

test_that("dem() truncates events and issues warning when exogenous_end is before last event", {
  # Create a small event set with events up to t=10
  events <- matrix(c(
    1, 1, 2, 1,
    2, 1, 2, 0,
    10, 2, 3, 1,
    11, 2, 3, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")
  
  # exogenous_end = 5 (before the last event)
  expect_warning(
    fit <- dem(events, n_nodes = 3, formula_0_1 = ~1, formula_1_0 = ~1, 
               exogenous_end = 5, control = control.redeem(return_data = TRUE)),
    "exogenous_end is before the last event. Truncating events."
  )
  
  # Check max_time is truncated
  expect_equal(fit$max_time, 5)
  
  # Check events in fit are truncated
  expect_true(all(fit$events[, 1] <= 5))
  expect_equal(nrow(fit$events), 3) # (1,1), (2,0) and the type-3 at 5
  
  # Check preprocessed data covers the expected time
  # model_0_1$data should have intervals ending at 5
  data_0_1 <- fit$model_0_1$data
  expect_equal(max(data_0_1$time_new), 5)
})

test_that("rem() truncates events and issues warning when exogenous_end is before last event", {
  events <- matrix(c(
    1, 1, 2, 1,
    10, 2, 3, 1
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")
  
  expect_warning(
    fit <- rem(events, n_nodes = 3, formula = ~1, 
               exogenous_end = 5, control = control.redeem(return_data = TRUE)),
    "exogenous_end is before the last event. Truncating events."
  )
  
  expect_equal(fit$max_time, 5)
  expect_true(all(fit$events[, 1] <= 5))
  
  data <- fit$model$data
  expect_equal(max(data$time_new), 5)
})

test_that("exogenous_end extends observation window if after last event", {
  events <- matrix(c(
    1, 1, 2, 1,
    2, 1, 2, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")
  
  # exogenous_end = 5 (after last event at 2)
  fit <- dem(events, n_nodes = 3, formula_0_1 = ~1, formula_1_0 = ~1, 
             exogenous_end = 5, control = control.redeem(return_data = TRUE))
  
  expect_equal(fit$max_time, 5)
  data_0_1 <- fit$model_0_1$data
  expect_equal(max(data_0_1$time_new), 5)
  
  # Check that there is an interval ending at 5
  expect_true(any(data_0_1$time_new == 5))
  expect_equal(max(data_0_1$time_new), 5)
})
