library(testthat)
library(redeem)

test_that("S3 methods work on actual model fits (not AI dummy mocks)", {
  # 1. Prepare small simulated dataset
  events <- matrix(c(
    1.0, 1, 2, 1,
    2.0, 1, 2, 0,
    3.0, 2, 3, 1,
    4.0, 2, 3, 0,
    5.0, 1, 3, 1,
    6.0, 1, 3, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")
  n_nodes <- 3

  # 2. Fit a simple DEM model (incidence and duration) via NR
  fit_dem_simple <- dem(
    events = events,
    formula_0_1 = ~ Intercept,
    formula_1_0 = ~ Intercept,
    n_nodes = n_nodes,
    directed = TRUE,
    control = control.redeem(estimate = "NR", save_hist = TRUE, check_matrix = FALSE)
  )

  # 3. Fit a complex DEM model via Blockwise (to estimate popularity & baseline changepoints)
  fit_dem_complex <- dem(
    events = events,
    formula_0_1 = ~ Intercept + degree(type = "out_sender") + baseline(changepoints = c(3.0)),
    formula_1_0 = ~ Intercept + degree(type = "out_receiver") + baseline(changepoints = c(3.0)),
    n_nodes = n_nodes,
    directed = TRUE,
    control = control.redeem(estimate = "Blockwise", save_hist = TRUE, it_max = 5, check_matrix = FALSE)
  )

  # 4. Fit a REM model
  fit_rem <- rem(
    events = events,
    formula = ~ Intercept + baseline(changepoints = c(3.0)),
    n_nodes = n_nodes,
    directed = TRUE,
    control = control.redeem(estimate = "NR", save_hist = TRUE, check_matrix = FALSE)
  )

  # 5. Fit a single-model DEM
  fit_dem_single <- dem(
    events = events,
    formula_0_1 = ~ Intercept,
    formula_1_0 = NULL,
    n_nodes = n_nodes,
    directed = TRUE,
    control = control.redeem(estimate = "NR", check_matrix = FALSE)
  )

  # --- Test Print & Summary ---
  expect_output(print(fit_dem_simple), "Durational Event Model")
  expect_output(print(fit_rem), "Relational Event Model")

  s_dem_simple <- summary(fit_dem_simple)
  expect_s3_class(s_dem_simple, "summary.dem")
  expect_output(print(s_dem_simple), "Results for Incidence Intensity")

  s_dem_complex <- summary(fit_dem_complex)
  expect_s3_class(s_dem_complex, "summary.dem")
  expect_output(print(s_dem_complex), "Degree Effects Summary|Degree Effects")

  s_dem_single <- summary(fit_dem_single)
  expect_s3_class(s_dem_single, "summary.dem")
  expect_true(is.na(s_dem_single$AIC))

  s_rem <- summary(fit_rem)
  expect_s3_class(s_rem, "summary.rem")
  expect_output(print(s_rem), "Log-likelihood:")

  # --- Test logLik ---
  expect_s3_class(logLik(fit_dem_simple$model_0_1), "logLik")
  expect_s3_class(logLik(fit_rem), "logLik")
  expect_equal(as.numeric(logLik(fit_rem)), fit_rem$model$llh)

  # --- Test Predict ---
  pred_dem <- predict(fit_dem_simple)
  expect_s3_class(pred_dem, "data.frame")
  expect_true(any(pred_dem$mode == "formation"))

  pred_rem <- predict(fit_rem)
  expect_s3_class(pred_rem, "data.frame")

  # --- Test Plotting ---
  pdf(NULL) # Disable plot output window
  on.exit(dev.off(), add = TRUE)

  # Plot models
  expect_silent(plot(fit_dem_simple))
  expect_silent(plot(fit_dem_complex))
  expect_silent(plot(fit_rem))

  # Plot baselines
  expect_silent(plot_baseline(fit_dem_complex, process = "formation"))
  expect_silent(plot_baseline(fit_dem_complex, process = "dissolution"))
  expect_silent(plot_baseline(fit_rem))

  # Plot residuals
  resids_dem <- get_residuals(fit_dem_simple)
  expect_silent(plot(resids_dem))
  
  resids_rem <- get_residuals(fit_rem)
  expect_silent(plot(resids_rem))
})

test_that("summary.redeem_result edge cases", {
  # 1. Prepare simple model fit
  events <- matrix(c(
    1.0, 1, 2, 1,
    2.0, 1, 2, 0,
    3.0, 2, 3, 1,
    4.0, 2, 3, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")
  
  fit <- dem(
    events = events,
    formula_0_1 = ~ Intercept,
    formula_1_0 = NULL,
    n_nodes = 3,
    directed = TRUE,
    control = control.redeem(estimate = "NR", check_matrix = FALSE)
  )
  
  model <- fit$model_0_1

  # Test negative variance-covariance
  model_neg <- model
  model_neg$covariance <- matrix(-1, 1, 1)
  s_neg <- summary(model_neg)
  expect_equal(s_neg$coefficients[1, "Std. Error"], 0)

  # Test NaN t-value
  model_nan <- model
  model_nan$est_core[1] <- 0
  model_nan$covariance <- matrix(0, 1, 1)
  s_nan <- summary(model_nan)
  expect_true(is.na(s_nan$coefficients[1, "t value"]))

  # Test no core effects
  model_none <- model
  model_none$est_core <- numeric(0)
  model_none$covariance <- NULL
  s_none <- summary(model_none)
  expect_output(print(s_none), "No fixed covariate effects.")

  # Test no history plotting warning
  fit_no_hist <- dem(
    events = events,
    formula_0_1 = ~ Intercept,
    formula_1_0 = NULL,
    n_nodes = 3,
    directed = TRUE,
    control = control.redeem(estimate = "NR", save_hist = FALSE, check_matrix = FALSE)
  )
  expect_message(plot(fit_no_hist$model_0_1), "Trace plots for coefficients are unavailable")
})

test_that("Ranking and GOF S3 methods coverage", {
  ranking <- data.frame(Cutpoint = 1:5, Recall = (1:5)/5, Precision = (5:1)/5)
  class(ranking) <- c("ranking_redeem", "data.frame")
  expect_output(print(ranking), "Ranking Results")
  
  # Long ranking
  ranking_long <- data.frame(Cutpoint = 1:15, Recall = (1:15)/15, Precision = (15:1)/15)
  class(ranking_long) <- c("ranking_redeem", "data.frame")
  expect_output(print(ranking_long), "truncated")

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  # Test plot metrics
  expect_silent(plot(ranking, metric = "recall"))
  expect_silent(plot(ranking, metric = "precision"))
  expect_silent(plot(ranking, metric = "both"))
})
