library(testthat)
library(redeem)

test_that("get_oos_likelihood() for REM matches the exact theoretical log-likelihood", {
    n_nodes <- 3

    train_events <- matrix(c(
        1.0, 1, 2, 1,
        2.0, 2, 3, 1
    ), ncol = 4, byrow = TRUE)
    colnames(train_events) <- c("time", "from", "to", "type")

    fit <- rem(events = train_events, formula = ~intercept, n_nodes = n_nodes, directed = TRUE)

    test_events <- matrix(c(
        10.0, 1, 3, 1
    ), ncol = 4, byrow = TRUE)

    ll <- get_oos_likelihood(fit, edgelist_test = test_events)

    # All 6 dyads are available, and all have the exact same LP.
    # Therefore, log-probability of any observed event is log(1/6)
    expect_equal(length(ll), 1)
    expect_equal(ll[1], log(1/6))
})

test_that("get_oos_likelihood() for DEM matches the exact theoretical log-likelihood under blocking", {
    n_nodes <- 3

    train_events <- matrix(c(
        1.0, 1, 2, 1,   # 1 and 2 start interacting
        2.0, 1, 2, 0,   # 1 and 2 end interacting
        3.0, 2, 3, 1,   # 2 and 3 start interacting
        4.0, 2, 3, 0    # 2 and 3 end interacting
    ), ncol = 4, byrow = TRUE)
    colnames(train_events) <- c("time", "from", "to", "type")

    # Fit a simple DEM model with intercepts
    fit <- dem(events = train_events, formula_0_1 = ~intercept, formula_1_0 = ~intercept, n_nodes = n_nodes, directed = TRUE, simultaneous_interactions = FALSE)

    # Test event is the formation at t=4.5 and termination at t=5.0
    test_events <- matrix(c(
        4.5, 1, 2, 1,
        5.0, 1, 2, 0
    ), ncol = 4, byrow = TRUE)

    ll <- get_oos_likelihood(fit, edgelist_test = test_events)

    # 1. At t=4.5: no active interactions. All 6 dyads are available for formation.
    # Probability is 1/6.
    # 2. At t=5.0: (1, 2) is active. Nodes 1 & 2 are busy, blocking all other dyads.
    # Probability of dissolution is 1.0.
    expect_equal(length(ll), 2)
    expect_equal(ll[1], log(1/6))
    expect_equal(ll[2], 0)
})

test_that("get_oos_likelihood() replaces -Inf degrees with the minimum of finite degrees", {
    n_nodes <- 4
    train_events <- matrix(c(
        1.0, 1, 2, 1,
        2.0, 2, 3, 1,
        3.0, 1, 2, 0,
        4.0, 2, 3, 0
    ), ncol = 4, byrow = TRUE)
    colnames(train_events) <- c("time", "from", "to", "type")

    # Fit a real REM model with intercept and degree
    fit <- rem(events = train_events, formula = ~intercept + degree, n_nodes = n_nodes, directed = FALSE)

    # Manually inject -Inf values to the estimated degree parameters
    fit$model$est_degree <- c(-Inf, -2.0, -1.5, -Inf)

    test_events <- matrix(c(
        5.0, 1, 2, 1
    ), ncol = 4, byrow = TRUE)

    # If -Inf was not replaced, the degree of node 1 (-Inf) would result in a -Inf log-likelihood.
    # By replacing -Inf with the minimum of finite degrees (-2.0), the log-likelihood stays finite.
    ll <- get_oos_likelihood(fit, edgelist_test = test_events)
    expect_true(all(is.finite(ll)))
})

test_that("predict_baseline_trend() correctly identifies daily and weekly seasonal patterns", {
    # 1. Test Daily seasonality (2 days of hourly data)
    times_d <- seq(from = 0, to = 2 * 86400, by = 3600)
    trend_d <- times_d / 86400 * 0.1
    cycle_d <- sin(2 * pi * times_d / 86400)
    y_d <- trend_d + cycle_d

    dummy_model_d <- list(
        est_time = y_d,
        time_changepoints = times_d[-1],
        full_baseline = TRUE
    )

    # Predict into day 3
    target_d <- seq(from = 2 * 86400 + 3600, to = 3 * 86400, by = 3600)
    pred_d <- predict_baseline_trend(dummy_model_d, target_times = target_d)

    expect_equal(length(pred_d), length(target_d))
    expect_true(all(is.finite(pred_d)))

    # 2. Test Weekly seasonality (8 days of hourly data)
    times_w <- seq(from = 0, to = 8 * 86400, by = 3600)
    trend_w <- times_w / 86400 * 0.05
    cycle_w <- sin(2 * pi * times_w / 604800) # weekly period
    y_w <- trend_w + cycle_w

    dummy_model_w <- list(
        est_time = y_w,
        time_changepoints = times_w[-1],
        full_baseline = TRUE
    )

    # Predict into day 9-10
    target_w <- seq(from = 8 * 86400 + 3600, to = 10 * 86400, by = 3600)
    pred_w <- predict_baseline_trend(dummy_model_w, target_times = target_w)

    expect_equal(length(pred_w), length(target_w))
    expect_true(all(is.finite(pred_w)))
})

test_that("get_oos_likelihood() supports different baseline methods", {
    n_nodes <- 3

    train_events <- matrix(c(
        1.0, 1, 2, 1,
        2.0, 2, 3, 1,
        3.0, 1, 3, 1
    ), ncol = 4, byrow = TRUE)
    colnames(train_events) <- c("time", "from", "to", "type")

    fit <- rem(events = train_events, formula = ~intercept, n_nodes = n_nodes, directed = TRUE)

    # Manually inject some baseline estimates
    fit$model$est_time <- c(0.1, 0.2, 0.3)
    fit$model$time_changepoints <- c(1.5, 2.5, 3.5)
    fit$model$full_baseline <- TRUE

    test_events <- matrix(c(
        4.0, 1, 3, 1
    ), ncol = 4, byrow = TRUE)

    # "last" should use 0.3
    ll_last <- get_oos_likelihood(fit, edgelist_test = test_events, baseline_method = "last")
    # "beginning" should use 0.0
    ll_beg <- get_oos_likelihood(fit, edgelist_test = test_events, baseline_method = "beginning")
    # "mean" should use mean(c(0.1, 0.2, 0.3)) = 0.2
    ll_mean <- get_oos_likelihood(fit, edgelist_test = test_events, baseline_method = "mean")
    # "trend" should use predict_baseline_trend extrapolation
    ll_trend <- suppressWarnings(get_oos_likelihood(fit, edgelist_test = test_events, baseline_method = "trend"))

    expect_true(is.finite(ll_last))
    expect_true(is.finite(ll_beg))
    expect_true(is.finite(ll_mean))
    expect_true(is.finite(ll_trend))
    
    # Check that they all equal the theoretical value log(1/6)
    expect_equal(ll_last[1], log(1/6))
    expect_equal(ll_beg[1], log(1/6))
    expect_equal(ll_mean[1], log(1/6))
    expect_equal(ll_trend[1], log(1/6))
})


