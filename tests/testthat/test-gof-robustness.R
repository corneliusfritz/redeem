library(testthat)
library(redeem)

test_that("get_residuals() is equivalent with and without return_data for REM", {
    set.seed(42)
    n_nodes <- 5

    # Simulate data
    events <- rem.simulate(
        formula = ~intercept,
        coef = -2,
        n_events = 50,
        n_nodes = n_nodes,
        directed = FALSE
    )

    # Fit with return_data = TRUE
    fit_with <- rem(
        events = events,
        formula = ~intercept,
        n_nodes = n_nodes,
        directed = FALSE,
        control = control.redeem(return_data = TRUE)
    )

    # Fit with return_data = FALSE
    fit_without <- rem(
        events = events,
        formula = ~intercept,
        n_nodes = n_nodes,
        directed = FALSE,
        control = control.redeem(return_data = FALSE)
    )

    # Check that fit_without doesn't have data
    expect_null(fit_without$model$data)
    expect_false(fit_without$return_data)

    # Get residuals (using raw = TRUE for mathematical equivalence check)
    resids_with <- get_residuals(fit_with, raw = TRUE)
    resids_without <- get_residuals(fit_without, raw = TRUE)

    # Compare
    expect_equal(resids_with, resids_without, tolerance = 1e-10)
})

test_that("get_residuals() is equivalent with and without return_data for DEM", {
    set.seed(42)
    n_nodes <- 5

    # Simulate data
    events <- dem.simulate(
        formula_0_1 = ~intercept,
        formula_1_0 = ~intercept,
        coef_0_1 = -1,
        coef_1_0 = -1.5,
        n_events = 50,
        n_nodes = n_nodes,
        directed = FALSE
    )

    # Close any open interactions at the end to satisfy dem() validation
    unique_dyads <- unique(events[, 2:3, drop=FALSE])
    for(i in seq_len(nrow(unique_dyads))) {
        d1 <- unique_dyads[i, 1]
        d2 <- unique_dyads[i, 2]
        dyad_events <- events[events[,2] == d1 & events[,3] == d2, , drop=FALSE]
        if(nrow(dyad_events) > 0 && dyad_events[nrow(dyad_events), 4] == 1) {
            events <- rbind(events, c(max(events[,1]) + 1, d1, d2, 0))
        }
    }
    events <- events[order(events[,1]), ]

    # Fit with return_data = TRUE
    fit_with <- dem(
        events = events,
        formula_0_1 = ~intercept,
        formula_1_0 = ~intercept,
        n_nodes = n_nodes,
        directed = FALSE,
        control = control.redeem(return_data = TRUE)
    )

    # Fit with return_data = FALSE
    fit_without <- dem(
        events = events,
        formula_0_1 = ~intercept,
        formula_1_0 = ~intercept,
        n_nodes = n_nodes,
        directed = FALSE,
        control = control.redeem(return_data = FALSE)
    )

    # Check data absence
    expect_null(fit_without$model_0_1$data)
    expect_null(fit_without$model_1_0$data)

    # Get residuals (using raw = TRUE for mathematical equivalence check)
    resids_with <- get_residuals(fit_with, raw = TRUE)
    resids_without <- get_residuals(fit_without, raw = TRUE)

    # Compare
    expect_equal(resids_with, resids_without, tolerance = 1e-10)
})

test_that("reproduce_model_data works with build_time for REM", {
    set.seed(42)
    n_nodes <- 5
    events <- rem.simulate(
        formula = ~intercept,
        coef = -2,
        n_events = 100,
        n_nodes = n_nodes,
        directed = FALSE
    )

    # Use build_time
    fit <- rem(
        events = events,
        formula = ~intercept,
        n_nodes = n_nodes,
        directed = FALSE,
        control = control.redeem(build_time = 0.5, return_data = FALSE)
    )

    expect_equal(fit$build_time, 0.5)

    # Should be able to get residuals
    resids <- get_residuals(fit, raw = TRUE)
    expect_true(nrow(resids$resid_0_1) > 0)

    # Verify that the reconstructed data respects build_time
    # Internal check:
    repro <- redeem:::reproduce_model_data(fit)
    expect_true(min(repro$time) >= 0.5)
})
