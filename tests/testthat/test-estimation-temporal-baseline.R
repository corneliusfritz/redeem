library(testthat)
library(redeem)

test_that("dem() handles model-specific changepoints via baseline()", {
  n <- 5
  # Simple events: dyad (1,2) starts at 1, ends at 10; starts at 20, ends at 30
  events <- matrix(c(
    1, 1, 2, 1,
    10, 1, 2, 0,
    20, 1, 2, 1,
    30, 1, 2, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")

  # Estimate models with different changepoints and labels
  fit <- dem(
    events = events,
    n_nodes = n,
    formula_0_1 = ~  baseline(changepoints = 15, labels = "Late"),
    formula_1_0 = ~ baseline(changepoints = c(5, 25), labels = c("Mid", "End")),
    control = control.redeem(estimate = "NR", it_max = 2)
  )

  expect_s3_class(fit, "dem")

  # Check if model coefficients include the expected temporal labels
  expect_true(any(grepl("Late", names(fit$model_0_1$coefficients))))
  expect_true(any(grepl("Mid", names(fit$model_1_0$coefficients))))
  expect_true(any(grepl("End", names(fit$model_1_0$coefficients))))
})

test_that("rem() handles changepoints via baseline()", {
  n <- 5
  events <- matrix(c(
    1, 1, 2, 1,
    20, 1, 2, 1
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")

  fit <- rem(
    events = events,
    n_nodes = n,
    formula = ~ baseline(changepoints = 15, labels = "After15"),
    control = control.redeem(estimate = "NR", it_max = 2)
  )

  expect_s3_class(fit, "rem")
  expect_true(any(grepl("After15", names(fit$model$coefficients))))
})

test_that("Blockwise temporal baseline updates guarantee strict monotonicity of log-likelihood", {
  set.seed(42)
  n_nodes <- 6
  events <- dem.simulate(
    formula_0_1 = ~ 1,
    formula_1_0 = ~ 1,
    coef_0_1 = -1.5,
    coef_1_0 = -2.0,
    n_events = 200,
    n_nodes = n_nodes,
    directed = TRUE
  )
  
  max_t <- max(events[, 1])
  changepoints <- seq(max_t / 10, max_t - max_t / 10, length.out = 9)
  
  fit <- dem(
    events = events,
    formula_0_1 = ~ baseline(changepoints = changepoints),
    formula_1_0 = ~ baseline(changepoints = changepoints),
    n_nodes = n_nodes,
    directed = TRUE,
    control = control.redeem(
      estimate = "Blockwise",
      tol = 1e-12,
      it_max = 20,
      save_hist = TRUE
    )
  )
  
  # Check 0-1 process likelihood history
  llh_hist_0_1 <- fit$model_0_1$llh_hist
  llh_hist_0_1 <- llh_hist_0_1[llh_hist_0_1 != 0 & !is.na(llh_hist_0_1)]
  if (length(llh_hist_0_1) > 1) {
    diffs_0_1 <- diff(llh_hist_0_1)
    expect_true(all(diffs_0_1 >= -1e-9), info = paste("LLH decreased in 0-1 process updates:", paste(diffs_0_1, collapse = ", ")))
  }
  
  # Check 1-0 process likelihood history
  llh_hist_1_0 <- fit$model_1_0$llh_hist
  llh_hist_1_0 <- llh_hist_1_0[llh_hist_1_0 != 0 & !is.na(llh_hist_1_0)]
  if (length(llh_hist_1_0) > 1) {
    diffs_1_0 <- diff(llh_hist_1_0)
    expect_true(all(diffs_1_0 >= -1e-9), info = paste("LLH decreased in 1-0 process updates:", paste(diffs_1_0, collapse = ", ")))
  }
})

