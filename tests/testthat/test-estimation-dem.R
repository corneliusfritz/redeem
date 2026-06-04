library(testthat)
library(redeem)

test_that("dem() works with simultaneous_interactions = TRUE", {
    set.seed(123)
    n_nodes <- 4
    events <- dem.simulate(
        formula_0_1 = ~intercept,
        formula_1_0 = ~intercept,
        coef_0_1 = -1,
        coef_1_0 = -1.5,
        n_events = 20,
        n_nodes = n_nodes,
        directed = FALSE
    )

    fit <- dem(
        events = events,
        formula_0_1 = ~intercept,
        formula_1_0 = ~intercept,
        n_nodes = n_nodes,
        directed = FALSE,
        simultaneous_interactions = TRUE
    )

    expect_s3_class(fit, "dem")
    expect_true(fit$simultaneous_interactions)
})

test_that("dem() handles empty results gracefully", {
    # Test with data that would generate very few transitions if possible
    events <- matrix(c(0, 1, 2, 0), ncol = 4)
    colnames(events) <- c("time", "from", "to", "event")

    # Should probably error or return a minimal object if no transitions are found
    expect_error(dem(events = events, formula_0_1 = ~1, formula_1_0 = ~1, n_nodes = 3))
})

test_that("dem() handles build_time argument through control", {
    set.seed(123)
    n_nodes <- 4
    events <- dem.simulate(
        formula_0_1 = ~intercept,
        formula_1_0 = ~intercept,
        coef_0_1 = -1,
        coef_1_0 = -1.5,
        n_events = 20,
        n_nodes = n_nodes,
        directed = FALSE
    )

    fit <- dem(
        events = events,
        formula_0_1 = ~intercept,
        formula_1_0 = ~intercept,
        n_nodes = n_nodes,
        directed = FALSE,
        control = control.redeem(build_time = 0.5)
    )

    expect_s3_class(fit, "dem")
})
