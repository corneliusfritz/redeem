library(testthat)
library(redeem)

test_that("get_ranking() works for windowed rem objects", {
    set.seed(123)
    n_nodes <- 5

    # Simulate data with a windowed effect
    # We use a simple intercept + windowed inertia
    events <- rem.simulate(
        formula = ~ intercept + inertia(window = 10),
        coef = c("intercept" = -3, "inertia" = 1),
        n_events = 100,
        n_nodes = n_nodes,
        directed = TRUE
    )

    # Split into train and test
    train_events <- events[1:80, ]
    test_events <- events[81:100, ]

    # Fit the model
    fit <- rem(
        events = train_events,
        formula = ~ intercept + inertia(window = 10),
        n_nodes = n_nodes,
        directed = TRUE
    )

    # Compute ranking
    ranking <- get_ranking(fit, edgelist_test = test_events, k_max = 10)

    expect_s3_class(ranking, "data.frame")
    expect_true(all(c("Cutpoint", "Recall") %in% names(ranking)))
    expect_equal(nrow(ranking), 11) # 0 to 10
    expect_true(all(ranking$Recall >= 0 & ranking$Recall <= 1))
})

test_that("get_ranking() works for windowed dem objects", {
    set.seed(123)
    n_nodes <- 5

    # Simulate data for dem with windowed effects
    events <- dem.simulate(
        formula_0_1 = ~ intercept + inertia(window = 10),
        formula_1_0 = ~intercept,
        coef_0_1 = c("intercept" = -3, "inertia" = 1),
        coef_1_0 = c("intercept" = -2),
        n_events = 100,
        n_nodes = n_nodes,
        directed = TRUE
    )

    # Split into train and test
    train_events <- events[1:80, ]
    test_events <- events[81:100, ]

    # Fit the model
    fit <- dem(
        events = train_events,
        formula_0_1 = ~ intercept + inertia(window = 10),
        formula_1_0 = ~intercept,
        n_nodes = n_nodes,
        directed = TRUE
    )

    ranking <- get_ranking(fit, edgelist_test = test_events, k_max = 10)

    expect_s3_class(ranking, "data.frame")
    expect_true(all(c("Cutpoint", "Recall") %in% names(ranking)))
    expect_equal(nrow(ranking), 11)
})

test_that("get_residuals() works for windowed rem objects", {
    set.seed(123)
    n_nodes <- 5

    events <- rem.simulate(
        formula = ~ intercept + inertia(window = 10),
        coef = c("intercept" = -3, "inertia" = 1),
        n_events = 50,
        n_nodes = n_nodes,
        directed = TRUE
    )

    # Fit the model
    fit <- rem(
        events = events,
        formula = ~ intercept + inertia(window = 10),
        n_nodes = n_nodes,
        directed = TRUE,
        control = control.redeem(return_data = TRUE)
    )
    resids <- get_residuals(fit)

    expect_s3_class(resids, "data.frame")
    expect_true(all(c("time", "surv", "lower", "upper", "theoretical") %in% names(resids)))
    expect_true(nrow(resids) > 0)
})

test_that("get_residuals() works for windowed dem objects", {
    set.seed(123)
    n_nodes <- 5

    events <- dem.simulate(
        formula_0_1 = ~ intercept + inertia(window = 10),
        formula_1_0 = ~intercept,
        coef_0_1 = c("intercept" = -3, "inertia" = 1),
        coef_1_0 = c("intercept" = -2),
        n_events = 50,
        n_nodes = n_nodes,
        directed = TRUE
    )

    # Fit the model
    fit <- dem(
        events = events,
        formula_0_1 = ~ intercept + inertia(window = 10),
        formula_1_0 = ~intercept,
        n_nodes = n_nodes,
        directed = TRUE,
        control = control.redeem(return_data = TRUE)
    )

    resids <- get_residuals(fit)

    expect_s3_class(resids, "data.frame")
    expect_true(all(c("time", "surv", "lower", "upper", "theoretical") %in% names(resids)))
    expect_true(nrow(resids) > 0)
})

test_that("get_ranking() correctly updates windowed history during prediction", {
    set.seed(123)
    n_nodes <- 3

    # We want to see if an event in the test set influences subsequent events in the test set.
    # Formula: ~ inertia(window = 10)
    # Event 1: (1, 2) at time 100
    # Event 2: (1, 2) at time 105
    # If history is updated, inertia for Event 2 will be 1.

    train_events <- matrix(c(
        10, 1, 3, 1 # dummy event to train on something
    ), ncol = 4, byrow = TRUE)
    colnames(train_events) <- c("time", "from", "to", "type")

    test_events <- matrix(c(
        100, 1, 2, 1,
        105, 1, 2, 1
    ), ncol = 4, byrow = TRUE)
    colnames(test_events) <- c("time", "from", "to", "type")

    # Get ranking using internal helper
    # We use a very high coefficient for inertia to ensure it becomes the top prediction
    ranking <- redeem:::get_probabilities_per_test_event(
        terms = "inertia_wt10",
        data_list = list(matrix(0, 0, 0)),
        transformations = "identity",
        n_nodes = n_nodes,
        verbose = FALSE,
        directed = TRUE,
        coef_0_1 = 10,
        coef_1_0 = 0,
        degree_coef_0_1 = rep(0, 2 * n_nodes),
        degree_coef_1_0 = rep(0, 2 * n_nodes),
        simultaneous_interactions = FALSE,
        edgelist_train = train_events,
        edgelist_test = test_events,
        k = 1,
        is_rem = TRUE,
        window_info = list("10" = 10)
    )

    # Event 2 at t=105
    pred_t105 <- ranking[[2]]$predicted
    obs_t105 <- ranking[[2]]$observed

    expect_equal(as.numeric(obs_t105[1, 1:2]), c(1, 2))
    # If updated, predicted[1, 1:2] should be (1, 2)
    # If NOT updated, predicted[1, 1:2] will be some other dyad (e.g. (1,3) which had an event in training)
    expect_equal(as.numeric(pred_t105[1, 1:2]), c(1, 2))
})
