library(testthat)
library(redeem)

test_that("Automated changepoint detection in rem.simulate", {
    n_nodes <- 20

    # Continuous covariate that switches at t=1.0
    cov1 <- matrix(0, n_nodes, n_nodes)
    cov2 <- matrix(1, n_nodes, n_nodes)

    continuous_cov <- list("0" = cov1, "1.0" = cov2)

    formula_rem <- ~ intercept + dyadic_cov(data = continuous_cov)
    coef_gt <- c(-1, 1) # Higher intercept

    set.seed(42)
    events <- rem.simulate(
        formula = formula_rem,
        coef = coef_gt,
        n_nodes = n_nodes,
        time = 2.0,
        verbose = FALSE
    )

    # Verify we got events in both intervals
    expect_true(any(events[, 1] < 1.0))
    expect_true(any(events[, 1] > 1.0))

    n_first <- sum(events[, 1] < 1.0)
    n_second <- sum(events[, 1] >= 1.0)

    expect_true(n_second > n_first)
})

test_that("Automated changepoint detection in dem.simulate", {
    n_nodes <- 20

    cov1 <- matrix(0, n_nodes, n_nodes)
    cov2 <- matrix(2, n_nodes, n_nodes)

    continuous_cov <- list("0" = cov1, "0.5" = cov2)

    set.seed(42)
    events <- dem.simulate(
        formula_0_1 = ~ intercept + dyadic_cov(data = continuous_cov),
        formula_1_0 = ~ intercept,
        coef_0_1 = c(-1, 0.5),
        coef_1_0 = c(-2),
        n_nodes = n_nodes,
        time = 1.0,
        verbose = FALSE
    )

    expect_true(any(events[, 1] < 0.5))
    expect_true(any(events[, 1] > 0.5))
})
