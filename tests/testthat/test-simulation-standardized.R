library(testthat)
library(redeem)

test_that("dem.simulate standardized baseline logic works", {
    set.seed(123)
    n_nodes <- 5
    
    # Case 1: ~Intercept (Implicitly NOT full_baseline)
    # User provides 1 coefficient for 1 changepoint. First slice is Intercept + 0.
    events <- dem.simulate(
        formula_0_1 = ~Intercept + baseline(changepoints = 10),
        formula_1_0 = ~Intercept + baseline(changepoints = 10),
        coef_0_1 = -1,
        coef_1_0 = -1.5,
        baseline_0_1 = 1, # Shift at t=10
        baseline_1_0 = 1, # Shift at t=10
        time = 20,
        n_nodes = n_nodes,
        verbose = TRUE
    )
    
    # We can't easily check internal states here without print_simulation_info
    # but we can check if it runs without error.
    expect_true(nrow(events) > 0)
    
    # Case 2: ~0 (Full baseline needed)
    # User provides 2 coefficients for 1 changepoint.
    events_full <- dem.simulate(
        formula_0_1 = ~0 + baseline(changepoints = 10),
        formula_1_0 = ~0 + baseline(changepoints = 10),
        baseline_0_1 = c(-1, 0), # First slice is -1, second is 0
        baseline_1_0 = c(-1.5, -0.5),
        time = 20,
        n_nodes = n_nodes,
        verbose = TRUE
    )
    expect_true(nrow(events_full) > 0)
    
    # Case 3: Both n_events and time specified should error
    expect_error(
        dem.simulate(
            formula_0_1 = ~Intercept,
            time = 10,
            n_events = 10,
            n_nodes = n_nodes
        ),
        "not both"
    )
})

test_that("Simulation past last changepoint works", {
    set.seed(123)
    n_nodes <- 5
    # CP at 10, simulate until 100
    events <- dem.simulate(
        formula_0_1 = ~Intercept + baseline(changepoints = 10),
        formula_1_0 = ~Intercept + baseline(changepoints = 10),
        coef_0_1 = -2, # Reasonable baseline
        coef_1_0 = -1.5,
        baseline_0_1 = 1, # Increase intensity after 10
        baseline_1_0 = 0,
        time = 100,
        n_nodes = n_nodes
    )
    expect_true(max(events[,1]) > 50)
})
