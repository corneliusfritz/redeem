library(testthat)
library(redeem)
library(survival)

test_that("rem() semiparametric path runs and returns a coxph object", {
    set.seed(123)
    n_nodes <- 5
    
    # Simulate some data using the parametric path
    events <- rem.simulate(
        formula = ~1,
        coef = -2,
        n_events = 100,
        n_nodes = n_nodes,
        directed = FALSE
    )

    # Estimate using the semiparametric path
    fit_semi <- rem(
        events = events,
        formula = ~intercept,
        n_nodes = n_nodes,
        directed = FALSE,
        semiparametric = TRUE
    )

    # Check that model is a coxph object
    expect_s3_class(fit_semi$model, "coxph")
    
    # Check that summary works
    s <- summary(fit_semi)
    expect_s3_class(s, "summary.rem")
    
    # Check that print doesn't crash
    expect_output(print(fit_semi), "Relational Event Model")
})

test_that("rem() semiparametric path handles covariates", {
    set.seed(42)
    n_nodes <- 6
    
    # Create dyadic covariate
    dyad_cov <- matrix(rnorm(n_nodes * n_nodes), n_nodes, n_nodes)
    dyad_cov <- (dyad_cov + t(dyad_cov)) / 2
    diag(dyad_cov) <- 0

    # Simulate data
    events <- rem.simulate(
        formula = ~1 + dyad_cov,
        coef = c(-1.5, 0.8),
        n_events = 200,
        n_nodes = n_nodes,
        directed = FALSE
    )

    # Estimate using the semiparametric path
    # Note: intercept is usually absorbed in baseline for Cox, but if explicitly requested 
    # as a covariate it might be handled differently depending on how formula_new is built.
    # In helper_estimation.R: update(formula_new, new = "... ~ . - Intercept")
    fit_semi <- rem(
        events = events,
        formula = ~ dyad_cov,
        n_nodes = n_nodes,
        directed = FALSE,
        semiparametric = TRUE
    )

    expect_s3_class(fit_semi$model, "coxph")
    # Coefficients should be proximal to true value (hazard ratios)
    expect_true(length(coef(fit_semi$model)) == 1)
})
