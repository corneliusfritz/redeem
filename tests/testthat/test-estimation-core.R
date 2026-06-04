library(testthat)
library(redeem)
library(data.table)

test_that("update_core_r, update_core_glm, and update_core_cpp are consistent", {
    set.seed(123)
    n <- 100
    p <- 3
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- paste0("x", 1:p)
    beta_true <- c(0.5, -0.2, 0.1)
    lambda <- exp(X %*% beta_true)
    y <- rpois(n, lambda)

    data <- data.table(event = y, pair_id = 1:n)
    covariates <- X
    est_core <- rep(0, p)
    identifiable <- rep(TRUE, p)
    offset_fixed <- rep(0, n)
    weights <- rep(1, n)

    # 1. Test update_core_r
    prediction <- exp(offset_fixed + covariates %*% est_core)
    est_r <- redeem:::update_core_r(
        data = data,
        covarites = covariates,
        prediction = prediction,
        est_core = est_core,
        identifiable = identifiable,
        offset_fixed = offset_fixed
    )

    # 2. Test update_core_glm
    est_glm <- redeem:::update_core_glm(
        data = data,
        covarites = covariates,
        est_core = est_core,
        identifiable = identifiable,
        offset_fixed = offset_fixed
    )

    # 3. Test update_core_cpp
    est_cpp <- redeem:::update_core_cpp(
        X = covariates,
        y = y,
        prediction = prediction,
        est_core = est_core,
        identifiable = which(identifiable) - 1,
        offset_fixed = offset_fixed,
        weights = weights
    )

    # Compare with standard GLM (convergence)
    fit_glm_full <- stats::glm(y ~ X - 1, family = poisson(), offset = offset_fixed)
    beta_glm <- as.numeric(coef(fit_glm_full))

    # Run update_core_r iteratively
    est_r_iter <- est_core
    for(i in 1:10) {
        pred_iter <- exp(offset_fixed + covariates %*% est_r_iter)
        est_r_iter <- redeem:::update_core_r(
            data = data,
            covarites = covariates,
            prediction = pred_iter,
            est_core = est_r_iter,
            identifiable = identifiable,
            offset_fixed = offset_fixed
        )
    }
    expect_equal(as.numeric(est_r_iter), beta_glm, tolerance = 1e-5)

    # Run update_core_cpp iteratively
    est_cpp_iter <- est_core
    for(i in 1:10) {
        pred_iter <- exp(offset_fixed + covariates %*% est_cpp_iter)
        est_cpp_iter <- redeem:::update_core_cpp(
            X = covariates,
            y = y,
            prediction = pred_iter,
            est_core = est_cpp_iter,
            identifiable = which(identifiable) - 1,
            offset_fixed = offset_fixed,
            weights = weights
        )
    }
    expect_equal(as.numeric(est_cpp_iter), beta_glm, tolerance = 1e-5)

    # update_core_glm should match glm.fit(maxit=1) exactly
    fit_glm_1 <- suppressWarnings(stats::glm.fit(x = X, y = y, family = poisson(), offset = offset_fixed, control = list(maxit = 1)))
    expect_equal(as.numeric(est_glm), as.numeric(fit_glm_1$coefficients), tolerance = 1e-7)

    # Check that all updates improved the likelihood from beta=0
    llh_init <- sum(dpois(y, exp(offset_fixed), log = TRUE))
    llh_r <- sum(dpois(y, exp(offset_fixed + covariates %*% est_r), log = TRUE))
    llh_glm <- sum(dpois(y, exp(offset_fixed + covariates %*% est_glm), log = TRUE))
    llh_cpp <- sum(dpois(y, exp(offset_fixed + covariates %*% est_cpp), log = TRUE))

    expect_gt(llh_r, llh_init)
    expect_gt(llh_glm, llh_init)
    expect_gt(llh_cpp, llh_init)
})

test_that("Accelerated estimation works and improves performance/convergence", {
    set.seed(123)
    # Use a small directed model where degree effects can be decoupled
    n_nodes <- 10
    events <- dem.simulate(
        formula_0_1 = ~ intercept,
        formula_1_0 = ~ intercept,
        coef_0_1 = c("intercept" = -3),
        coef_1_0 = c("intercept" = -2),
        n_events = 200,
        n_nodes = n_nodes,
        directed = TRUE
    )

    # Standard estimation
    fit_std <- dem(
        events = events,
        formula_0_1 = ~ intercept,
        formula_1_0 = ~ intercept,
        n_nodes = n_nodes,
        directed = TRUE,
        control = control.redeem(accelerated = FALSE, it_max = 20)
    )

    # Accelerated estimation
    fit_acc <- dem(
        events = events,
        formula_0_1 = ~ intercept,
        formula_1_0 = ~ intercept,
        n_nodes = n_nodes,
        directed = TRUE,
        control = control.redeem(accelerated = TRUE, it_max = 20)
    )

    # Check coefficients are similar
    expect_equal(fit_std$model_0_1$coefficients, fit_acc$model_0_1$coefficients, tolerance = 1e-3)
    expect_equal(fit_std$model_0_1$est_degree, fit_acc$model_0_1$est_degree, tolerance = 1e-3)

    # In some cases acceleration might reach a better LLH or reach it faster
    expect_true(fit_acc$model_0_1$llh >= fit_std$model_0_1$llh - 1e-7)
})

test_that("Accelerated estimation works for undirected models", {
    set.seed(123)
    n_nodes <- 10
    # Simulate undirected data
    events <- dem.simulate(
        formula_0_1 = ~ intercept,
        formula_1_0 = ~ intercept,
        coef_0_1 = c("intercept" = -3),
        coef_1_0 = c("intercept" = -2),
        n_events = 200,
        n_nodes = n_nodes,
        directed = FALSE
    )

    # Standard estimation
    fit_std <- dem(
        events = events,
        formula_0_1 = ~ intercept,
        formula_1_0 = ~ intercept,
        n_nodes = n_nodes,
        directed = FALSE,
        control = control.redeem(accelerated = FALSE, it_max = 20)
    )

    # Accelerated estimation
    fit_acc <- dem(
        events = events,
        formula_0_1 = ~ intercept,
        formula_1_0 = ~ intercept,
        n_nodes = n_nodes,
        directed = FALSE,
        control = control.redeem(accelerated = TRUE, it_max = 20)
    )

    # Check coefficients are similar
    expect_equal(fit_std$model_0_1$coefficients, fit_acc$model_0_1$coefficients, tolerance = 1e-3)
    expect_equal(fit_std$model_0_1$est_degree, fit_acc$model_0_1$est_degree, tolerance = 1e-3)

    expect_true(fit_acc$model_0_1$llh >= fit_std$model_0_1$llh - 1e-7)
})

test_that("Accelerated estimation with time-varying baseline (estimate_mmt) works", {
    set.seed(123)
    n_nodes <- 8
    # Define a simple time-varying baseline
    events <- dem.simulate(
        formula_0_1 = ~ intercept + baseline(changepoints = 50),
        formula_1_0 = ~ intercept,
        coef_0_1 = c("intercept" = -3, "baseline_cat50" = 0.5),
        coef_1_0 = c("intercept" = -2),
        n_events = 150,
        n_nodes = n_nodes,
        directed = TRUE
    )

    # Standard estimation (estimate_mmt)
    fit_std <- dem(
        events = events,
        formula_0_1 = ~ intercept + baseline(changepoints = 50),
        formula_1_0 = ~ intercept,
        n_nodes = n_nodes,
        directed = TRUE,exogenous_end = 100,
        control = control.redeem(accelerated = FALSE, it_max = 20)
    )

    # Accelerated estimation (estimate_mmt)
    fit_acc <- dem(
        events = events,
        formula_0_1 = ~ intercept + baseline(changepoints = 50),
        formula_1_0 = ~ intercept,
        n_nodes = n_nodes,
        directed = TRUE,exogenous_end = 100,
        control = control.redeem(accelerated = TRUE, it_max = 20)
    )

    expect_equal(fit_std$model_0_1$coefficients, fit_acc$model_0_1$coefficients, tolerance = 1e-3)
    expect_true(fit_acc$model_0_1$llh >= fit_std$model_0_1$llh - 1e-7)
})

test_that("update_core_glm fallback retains correct dimensions with non-identifiable parameters", {
    set.seed(123)
    n <- 100
    p <- 3
    X <- matrix(rnorm(n * p), n, p)
    X[, 2] <- 1.0 # non-identifiable
    colnames(X) <- paste0("x", 1:p)
    beta_true <- c(0.5, 0.0, 0.1)
    lambda <- exp(X %*% beta_true)
    y <- rpois(n, lambda)

    data <- data.table(event = y, pair_id = 1:n)
    covariates <- X
    est_core <- beta_true
    identifiable <- c(TRUE, FALSE, TRUE)
    offset_fixed <- rep(0, n)

    # Force fallback block trigger with subsample and check dimension safety
    est_glm <- redeem:::update_core_glm(
        data = data,
        covarites = covariates,
        est_core = est_core,
        identifiable = identifiable,
        offset_fixed = offset_fixed,
        subsample = 0.1
    )

    expect_length(est_glm, p)
    expect_equal(est_glm[2], 0)
})

