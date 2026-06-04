library(testthat)
library(redeem)

test_that("S3 predict methods work for rem objects", {
  set.seed(123)
  n_nodes <- 5
  events <- rem.simulate(formula = ~1, coef = -2, n_events = 50, n_nodes = n_nodes)
  fit <- rem(events = events, formula = ~intercept, n_nodes = n_nodes, control = control.redeem(return_data = TRUE))

  # Test basic predict types
  pred_response <- predict(fit, type = "response")
  pred_lp <- predict(fit, type = "lp")
  pred_terms <- predict(fit, type = "terms")

  expect_s3_class(pred_response, "data.frame")
  expect_s3_class(pred_lp, "data.frame")
  expect_s3_class(pred_terms, "data.frame")

  expect_equal(names(pred_response), c("from", "to", "prediction", "mode"))
  expect_equal(names(pred_lp), c("from", "to", "prediction", "mode"))
  expect_equal(names(pred_terms), c("from", "to", "mode", "Intercept"))

  expect_equal(nrow(pred_response), nrow(fit$model$data))
  expect_equal(nrow(pred_lp), nrow(fit$model$data))
  expect_equal(nrow(pred_terms), nrow(fit$model$data))

  # Test mathematical consistency: response = exp(lp)
  expect_equal(pred_response$prediction, exp(pred_lp$prediction), tolerance = 1e-7)

  # Test filtering by a specific time point
  pred_time <- predict(fit, time = 1.0, type = "response")
  expect_s3_class(pred_time, "data.frame")
  expect_true(nrow(pred_time) <= nrow(pred_response))
  expect_equal(names(pred_time), c("from", "to", "prediction", "mode"))

  # Verify for type = "terms" at a specific time point
  pred_time_terms <- predict(fit, time = 1.0, type = "terms")
  expect_s3_class(pred_time_terms, "data.frame")
  expect_equal(names(pred_time_terms), c("from", "to", "mode", "Intercept"))

  # Test on-the-fly data reconstruction (when return_data = FALSE)
  fit_no_data <- rem(events = events, formula = ~intercept, n_nodes = n_nodes, control = control.redeem(return_data = FALSE))
  expect_null(fit_no_data$model$data)

  pred_reconstructed <- predict(fit_no_data, type = "response")
  expect_s3_class(pred_reconstructed, "data.frame")
  expect_equal(nrow(pred_reconstructed), nrow(pred_response))
})

test_that("S3 predict methods work for dem objects", {
  set.seed(123)
  n_nodes <- 5
  events <- dem.simulate(
    formula_0_1 = ~intercept,
    formula_1_0 = ~intercept,
    coef_0_1 = -1,
    coef_1_0 = -1.5,
    n_events = 50,
    n_nodes = n_nodes
  )
  
  fit <- dem(
    events = events,
    formula_0_1 = ~intercept,
    formula_1_0 = ~intercept,
    n_nodes = n_nodes,
    control = control.redeem(return_data = TRUE)
  )

  # Test process selection and types
  pred_both <- predict(fit, type = "response", process = "both")
  expect_s3_class(pred_both, "data.frame")
  expect_equal(names(pred_both), c("from", "to", "prediction", "mode"))

  pred_form <- predict(fit, type = "response", process = "formation")
  pred_diss <- predict(fit, type = "response", process = "dissolution")

  expect_s3_class(pred_form, "data.frame")
  expect_s3_class(pred_diss, "data.frame")

  # Combined output should be row-bound form + diss
  expect_equal(pred_both, rbind(pred_form, pred_diss))

  # Test terms
  pred_terms_form <- predict(fit, type = "terms", process = "formation")
  expect_s3_class(pred_terms_form, "data.frame")
  expect_equal(names(pred_terms_form), c("from", "to", "mode", "Intercept"))

  # Test time filtering
  pred_time_form <- predict(fit, time = 2.0, type = "response", process = "formation")
  expect_s3_class(pred_time_form, "data.frame")
  expect_equal(names(pred_time_form), c("from", "to", "prediction", "mode"))

  pred_time_terms_form <- predict(fit, time = 2.0, type = "terms", process = "formation")
  expect_s3_class(pred_time_terms_form, "data.frame")
  expect_equal(names(pred_time_terms_form), c("from", "to", "mode", "Intercept"))

  # Test combined time filtering (process = "both")
  pred_time_both <- predict(fit, time = 2.0, type = "response", process = "both")
  expect_s3_class(pred_time_both, "data.frame")
  expect_equal(names(pred_time_both), c("from", "to", "prediction", "mode"))
  
  # Check that it contains dyads from both formation and dissolution
  pred_time_form_only <- predict(fit, time = 2.0, type = "response", process = "formation")
  pred_time_diss_only <- predict(fit, time = 2.0, type = "response", process = "dissolution")
  expect_equal(nrow(pred_time_both), nrow(pred_time_form_only) + nrow(pred_time_diss_only))
  
  pred_time_terms_both <- predict(fit, time = 2.0, type = "terms", process = "both")
  expect_s3_class(pred_time_terms_both, "data.frame")
  expect_true("formation_Intercept" %in% names(pred_time_terms_both))
  expect_true("dissolution_Intercept" %in% names(pred_time_terms_both))

  # Test on-the-fly reconstruction
  fit_no_data <- dem(
    events = events,
    formula_0_1 = ~intercept,
    formula_1_0 = ~intercept,
    n_nodes = n_nodes,
    control = control.redeem(return_data = FALSE)
  )
  
  expect_null(fit_no_data$model_0_1$data)
  expect_null(fit_no_data$model_1_0$data)

  pred_reconstructed <- predict(fit_no_data, type = "response", process = "both")
  expect_s3_class(pred_reconstructed, "data.frame")
  expect_equal(nrow(pred_reconstructed), nrow(pred_both))
})

test_that("S3 predict methods are robust to missing metadata and case-sensitivity", {
  set.seed(123)
  n_nodes <- 5
  events <- rem.simulate(formula = ~1, coef = -2, n_events = 20, n_nodes = n_nodes)
  fit <- rem(events = events, formula = ~intercept, n_nodes = n_nodes, control = control.redeem(return_data = TRUE))

  # 1. Empty prediction rows returns empty data.frame with correct columns/length
  pred_empty_terms <- predict(fit, time = 9999.0, type = "terms")
  expect_s3_class(pred_empty_terms, "data.frame")
  expect_equal(nrow(pred_empty_terms), 0)
  expect_equal(names(pred_empty_terms), c("from", "to", "mode", "Intercept"))

  pred_empty_lp <- predict(fit, time = 9999.0, type = "lp")
  expect_s3_class(pred_empty_lp, "data.frame")
  expect_equal(nrow(pred_empty_lp), 0)
  expect_equal(names(pred_empty_lp), c("from", "to", "prediction", "mode"))

  # 2. Robust to missing model$directed and model$n_nodes
  fit_corrupted <- fit
  fit_corrupted$model$directed <- NULL
  fit_corrupted$model$n_nodes <- NULL
  
  pred_robust <- predict(fit_corrupted, type = "response")
  expect_equal(nrow(pred_robust), nrow(fit$model$data))

  # 3. Robust to intercept naming (e.g. data having lowercase 'intercept')
  corrupted_data <- fit$model$data
  names(corrupted_data)[names(corrupted_data) == "Intercept"] <- "intercept"
  fit_corrupted_data <- fit
  fit_corrupted_data$model$data <- corrupted_data
  
  pred_case <- predict(fit_corrupted_data, type = "response")
  expect_equal(pred_case, predict(fit, type = "response"))

  # 4. Robust to out-of-bounds / NA node indices (should fallback to 0 for degree effects)
  fit_deg <- fit
  fit_deg$model$est_degree <- c(0.1, 0.2, 0.3, 0.4, 0.5) # length 5
  fit_deg$model$directed <- FALSE
  fit_deg$model$n_nodes <- 5

  data_oob <- copy(fit$model$data)
  # Set some node indices out-of-bounds (e.g. 6) and NA
  data_oob[1, from := 6]
  data_oob[2, to := NA]

  pred_oob <- predict(fit_deg, type = "lp")
  expect_true(all(!is.na(pred_oob$prediction)))
})

test_that("S3 predict sums up time effects in terms option", {
  set.seed(123)
  n_nodes <- 5
  events <- rem.simulate(formula = ~1, coef = -2, n_events = 50, n_nodes = n_nodes)
  
  # Fit a model
  fit <- rem(events = events, formula = ~intercept, n_nodes = n_nodes, control = control.redeem(return_data = TRUE))
  
  # Mock time effects in the fit object
  fit_mock <- fit
  fit_mock$model$coefficients <- c("Intercept" = -2.0, "time_1" = 0.5)
  
  # Predict terms
  pred_terms <- predict(fit_mock, type = "terms")
  expect_s3_class(pred_terms, "data.frame")
  expect_true("time_effects" %in% names(pred_terms))
  expect_false("time_1" %in% names(pred_terms))
})

test_that("S3 predictions mathematically agree with each other and data prediction column", {
  set.seed(123)
  n_nodes <- 5
  events <- rem.simulate(formula = ~1, coef = -2, n_events = 50, n_nodes = n_nodes)
  
  # 1. Parametric Model
  fit_p <- rem(events = events, formula = ~intercept, n_nodes = n_nodes, control = control.redeem(return_data = TRUE))
  pred_resp_p <- predict(fit_p, type = "response")
  pred_lp_p <- predict(fit_p, type = "lp")
  pred_terms_p <- predict(fit_p, type = "terms")
  
  expect_equal(pred_resp_p$prediction, exp(pred_lp_p$prediction), tolerance = 1e-7)
  term_cols_p <- setdiff(names(pred_terms_p), c("from", "to", "mode"))
  expect_equal(unname(rowSums(pred_terms_p[, term_cols_p, drop = FALSE])), unname(pred_lp_p$prediction), tolerance = 1e-7)
  
  # data$prediction vs response * duration
  duration_p <- exp(fit_p$model$data$offset)
  expect_equal(fit_p$model$data$prediction, pred_resp_p$prediction * duration_p, tolerance = 1e-7)
  
  # 2. Blockwise Model
  fit_b <- rem(events = events, formula = ~baseline(changepoints = c(1.0, 2.0, 3.0)), n_nodes = n_nodes, control = control.redeem(return_data = TRUE))
  pred_resp_b <- predict(fit_b, type = "response")
  pred_lp_b <- predict(fit_b, type = "lp")
  pred_terms_b <- predict(fit_b, type = "terms")
  
  expect_equal(pred_resp_b$prediction, exp(pred_lp_b$prediction), tolerance = 1e-7)
  term_cols_b <- setdiff(names(pred_terms_b), c("from", "to", "mode"))
  expect_equal(unname(rowSums(pred_terms_b[, term_cols_b, drop = FALSE])), unname(pred_lp_b$prediction), tolerance = 1e-7)
  
  # Reconstruct integrated intensity using baseline integration
  data_b <- fit_b$model$data
  est_time <- fit_b$model$est_time
  time_changepoints <- fit_b$model$time_changepoints
  
  time_slices_from <- findInterval(data_b$time, c(-Inf, time_changepoints, Inf))
  time_slices_to <- findInterval(data_b$time_new, c(-Inf, time_changepoints, Inf))
  
  baseline_mult <- as.vector(redeem:::get_time_offset(
    from_slice_r = time_slices_from,
    to_slice_r = time_slices_to,
    from_time_r = data_b$time,
    to_time_r = data_b$time_new,
    est_time_r = if (isTRUE(fit_b$model$full_baseline)) exp(as.vector(est_time)) else exp(c(0, as.vector(est_time))),
    changepoints_r = c(time_changepoints, max(data_b$time_new, na.rm = TRUE) + 1e-6)
  ))
  
  expect_equal(fit_b$model$data$prediction, baseline_mult, tolerance = 1e-7)
})

