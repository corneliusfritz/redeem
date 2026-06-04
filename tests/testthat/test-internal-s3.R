library(testthat)
library(redeem)

test_that("S3 methods work for rem objects", {
    set.seed(123)
    n_nodes <- 5
    events <- rem.simulate(formula = ~1, coef = -2, n_events = 50, n_nodes = n_nodes)
    fit <- rem(events = events, formula = ~intercept, n_nodes = n_nodes)

    # debugonce(summary)
    s <- summary(fit)
    expect_s3_class(s, "summary.rem")
    expect_true("coefficients" %in% names(s))

    # logLik
    ll <- logLik(fit)
    expect_s3_class(ll, "logLik")
    expect_true(is.numeric(ll))

    # plot (check if it doesn't error)
    # expect_error(plot(fit), NA) # skip for now to avoid potential plot device issues
})

test_that("S3 methods work for dem objects", {
    set.seed(123)
    n_nodes <- 5
    events <- dem.simulate(formula_0_1 = ~intercept, formula_1_0 = ~intercept, coef_0_1 = -1, coef_1_0 = -1.5, n_events = 50, n_nodes = n_nodes)
    fit <- dem(events = events, formula_0_1 = ~intercept, formula_1_0 = ~intercept, n_nodes = n_nodes)

    # summary
    s <- summary(fit)
    expect_s3_class(s, "summary.dem")
    expect_true("AIC" %in% names(s))
    expect_true("BIC" %in% names(s))

    # logLik
    ll_0_1 <- logLik(fit$model_0_1)
    expect_s3_class(ll_0_1, "logLik")

    # plot
    # expect_error(plot(fit), NA)
})

test_that("S3 methods work for dem.mm objects", {

    set.seed(123)
    n_nodes <- 5
    true_coef_0_1 <- -1
    true_coef_1_0 <- -1.5

    events <- dem.simulate(
        formula_0_1 = ~Intercept,
        formula_1_0 = ~Intercept,
        coef_0_1 = true_coef_0_1,
        coef_1_0 = true_coef_1_0,
        n_events = 100,
        n_nodes = n_nodes
    )
    fit <- dem(events = events, formula_0_1 = ~Intercept, formula_1_0 = ~Intercept, n_nodes = n_nodes)

    expect_s3_class(fit$model_0_1, "dem.nr")

    # summary
    s <- summary(fit$model_0_1)
    expect_s3_class(s, "summary.redeem_result")

    # logLik
    ll <- logLik(fit$model_0_1)
    expect_s3_class(ll, "logLik")
})
