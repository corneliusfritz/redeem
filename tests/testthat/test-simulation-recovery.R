library(testthat)
library(redeem)

test_that("rem() parameter recovery: Undirected (intercept only)", {
    set.seed(123)
    n_nodes <- 5
    true_intercept <- -2

    # Simulate data
    events <- rem.simulate(
        formula = ~1,
        coef = true_intercept,
        n_events = 500,
        n_nodes = n_nodes,
        directed = FALSE
    )

    # Estimate model
    fit <- rem(
        events = events,
        formula = ~intercept,
        n_nodes = n_nodes,
        directed = FALSE,
        control = control.redeem(weighting = FALSE)
    )

    # Check recovery
    expect_equal(as.numeric(fit$model$coefficients[1]), true_intercept, tolerance = 0.3)
    expect_false(fit$directed)
})

test_that("rem() parameter recovery: Directed (intercept only)", {
    set.seed(123)
    n_nodes <- 5
    true_intercept <- -2

    # Simulate data
    events <- rem.simulate(
        formula = ~intercept,
        coef = true_intercept,
        n_events = 500,
        n_nodes = n_nodes,
        directed = TRUE
    )

    # Estimate model
    fit <- rem(
        events = events,
        formula = ~intercept,
        n_nodes = n_nodes,
        directed = TRUE
    )

    expect_equal(as.numeric(fit$model$coefficients[1]), true_intercept, tolerance = 0.3)
    expect_true(fit$directed)
})

test_that("dem() parameter recovery: Undirected (intercept only)", {
    set.seed(123)
    n_nodes <- 5
    true_coef_0_1 <- -1
    true_coef_1_0 <- -1.5

    # Simulate data
    events <- dem.simulate(
        formula_0_1 = ~1,
        formula_1_0 = ~1,
        coef_0_1 = true_coef_0_1,
        coef_1_0 = true_coef_1_0,
        n_events = 500,
        n_nodes = n_nodes,
        directed = FALSE
    )

    # Estimate model
    fit <- dem(
        events = events,
        formula_0_1 = ~intercept,
        formula_1_0 = ~intercept,
        n_nodes = n_nodes,
        directed = FALSE,
        simultaneous_interactions = FALSE
    )

    expect_equal(as.numeric(fit$model_0_1$coefficients[1]), true_coef_0_1, tolerance = 0.3)
    expect_equal(as.numeric(fit$model_1_0$coefficients[1]), true_coef_1_0, tolerance = 0.3)
    expect_false(fit$directed)
})

test_that("dem() parameter recovery: Directed (intercept only)", {
    set.seed(123)
    n_nodes <- 5
    true_coef_0_1 <- -1
    true_coef_1_0 <- -1.5

    # Simulate data
    events <- dem.simulate(
        formula_0_1 = ~intercept,
        formula_1_0 = ~intercept,
        coef_0_1 = true_coef_0_1,
        coef_1_0 = true_coef_1_0,
        n_events = 500,
        n_nodes = n_nodes,
        directed = TRUE
    )

    # Estimate model
    fit <- dem(
        events = events,
        formula_0_1 = ~intercept,
        formula_1_0 = ~intercept,
        n_nodes = n_nodes,
        directed = TRUE,
        simultaneous_interactions = FALSE
    )

    expect_equal(as.numeric(fit$model_0_1$coefficients[1]), true_coef_0_1, tolerance = 0.3)
    expect_equal(as.numeric(fit$model_1_0$coefficients[1]), true_coef_1_0, tolerance = 0.3)
    expect_true(fit$directed)
})

test_that("rem() parameter recovery: Undirected with covariate", {
    set.seed(123)
    n_nodes <- 5
    # Intercept and one covariate
    true_coefs <- c(-2, 0.5)

    # Create a symmetric dyadic covariate
    dyad_cov <- matrix(rnorm(n_nodes * n_nodes), n_nodes, n_nodes)
    dyad_cov <- (dyad_cov + t(dyad_cov)) / 2
    diag(dyad_cov) <- 0

    # Simulate data
    events <- rem.simulate(
        formula = ~1 + dyad_cov,
        coef = true_coefs,
        n_events = 1000,
        n_nodes = n_nodes,
        directed = FALSE
    )

    # Estimate model
    fit <- rem(
        events = events,
        formula = ~ intercept + dyad_cov,
        n_nodes = n_nodes,
        directed = FALSE
    )

    expect_equal(as.numeric(fit$model$coefficients), true_coefs, tolerance = 0.3)
})

test_that("rem() parameter recovery: Directed with covariate", {
    set.seed(123)
    n_nodes <- 5
    true_coefs <- c(-2, 0.5)

    # Create a non-symmetric dyadic covariate
    dyad_cov <- matrix(rnorm(n_nodes * n_nodes), n_nodes, n_nodes)
    diag(dyad_cov) <- 0

    # Simulate data
    events <- rem.simulate(
        formula = ~ 1 + dyad_cov,
        coef = true_coefs,
        n_events = 1000,
        n_nodes = n_nodes,
        directed = TRUE
    )

    # Estimate model
    fit <- rem(
        events = events,
        formula = ~ intercept + dyad_cov,
        n_nodes = n_nodes,
        directed = TRUE
    )

    expect_equal(as.numeric(fit$model$coefficients), true_coefs, tolerance = 0.3)
})

test_that("dem() parameter recovery: Undirected with covariate", {
    set.seed(123)
    n_nodes <- 5
    true_coef_0_1 <- c(-1, 0.5)
    true_coef_1_0 <- c(-1.5, -0.3)

    dyad_cov <- matrix(rnorm(n_nodes * n_nodes), n_nodes, n_nodes)
    dyad_cov <- (dyad_cov + t(dyad_cov)) / 2
    diag(dyad_cov) <- 0

    # Simulate data
    events <- dem.simulate(
        formula_0_1 = ~intercept + dyad_cov,
        formula_1_0 = ~intercept + dyad_cov,
        coef_0_1 = true_coef_0_1,
        coef_1_0 = true_coef_1_0,
        n_events = 1000,
        n_nodes = n_nodes,
        directed = FALSE
    )

    # Estimate model
    fit <- dem(
        events = events,
        formula_0_1 = ~intercept+ dyad_cov,
        formula_1_0 = ~intercept + dyad_cov,
        n_nodes = n_nodes,
        directed = FALSE,
        simultaneous_interactions = FALSE
    )

    expect_equal(as.numeric(fit$model_0_1$coefficients), true_coef_0_1, tolerance = 0.3)
    expect_equal(as.numeric(fit$model_1_0$coefficients), true_coef_1_0, tolerance = 0.3)
})

test_that("dem() parameter recovery: Directed with covariate", {
    set.seed(123)
    n_nodes <- 10
    true_coef_0_1 <- c(-1, 0.5)
    true_coef_1_0 <- c(-1.5, -0.3)

    dyad_cov <- matrix(rnorm(n_nodes * n_nodes), n_nodes, n_nodes)
    diag(dyad_cov) <- 0

    # Simulate data
    events <- dem.simulate(
        formula_0_1 = ~intercept + dyad_cov,
        formula_1_0 = ~intercept + dyad_cov,
        coef_0_1 = true_coef_0_1,
        coef_1_0 = true_coef_1_0,
        n_events = 3000,
        n_nodes = n_nodes,
        directed = TRUE
    )

    # Estimate model
    fit <- dem(
        events = events,
        formula_0_1 = ~ intercept + dyad_cov,
        formula_1_0 = ~ intercept + dyad_cov,
        n_nodes = n_nodes,
        directed = TRUE,
        simultaneous_interactions = FALSE
    )

    expect_equal(as.numeric(fit$model_0_1$coefficients), true_coef_0_1, tolerance = 0.3)
    expect_equal(as.numeric(fit$model_1_0$coefficients), true_coef_1_0, tolerance = 0.3)
})

test_that("dem() parameter recovery: Directed with degree", {
    set.seed(42)
    n_nodes_test <- 10

    degree_out <- rnorm(n_nodes_test, 0, 0.1)
    degree_in <- rnorm(n_nodes_test, 0, 0.1)
    degree_in[1] <- 0
    degree_0_1 <- c(degree_out, degree_in)

    formula_0_1 <- ~ degrees
    formula_1_0 <- ~ degrees
    # debugonce(dem.simulate)
    sim_data <- dem.simulate(
      formula_0_1 = formula_0_1,
      formula_1_0 = formula_1_0,
      coef_degree_0_1 =  degree_0_1,
      coef_degree_1_0 = degree_0_1,
      n_nodes = n_nodes_test,
      n_events = 5000,
      directed = TRUE, verbose = FALSE
    )
    fit <- dem(events = sim_data, formula_0_1 = formula_0_1, formula_1_0 = formula_1_0,
               n_nodes = n_nodes_test, directed = TRUE,
               control = control.redeem(verbose = FALSE))


    cor_sender_0_1 <- cor(degree_out, fit$model_0_1$est_degree[1:n_nodes_test])
    cor_receiver_0_1 <- cor(degree_in, fit$model_0_1$est_degree[(n_nodes_test+1):(2*n_nodes_test)])
    expect_true(cor_sender_0_1 > 0.75)
    expect_true(cor_receiver_0_1 > 0.75)
})
