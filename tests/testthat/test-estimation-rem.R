library(testthat)
library(redeem)

test_that("rem() function works with directed = FALSE (default)", {
    # Simple test case with undirected events
    # time, from, to, event
    events <- matrix(c(
        1, 1, 2, 1,
        2, 2, 3, 1,
        3, 3, 1, 1
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "event")

    # Basic rem call
    # Assuming 4 nodes and a simple intercept model
    fit <- rem(
        events = events,
        formula = ~intercept,
        n_nodes = 4,
        directed = FALSE
    )

    expect_s3_class(fit, "rem")
    expect_false(fit$directed)
})

test_that("rem() function works with directed = TRUE", {
    events <- matrix(c(
        1, 1, 2, 1,
        2, 2, 3, 1,
        3, 3, 1, 1
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "event")

    fit <- rem(
        events = events,
        formula = ~intercept,
        n_nodes = 4,
        directed = TRUE
    )

    expect_s3_class(fit, "rem")
    expect_true(fit$directed)
})

test_that("rem() handles missing n_nodes gracefully by automatically identifying it", {
    events <- matrix(c(1, 1, 2, 1), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "event")

    # n_nodes is now automatically identified
    fit <- rem(events = events, formula = ~1, directed = TRUE)
    expect_s3_class(fit, "rem")
    expect_equal(fit$n_nodes, 2)
})

test_that("rem() works with a single event", {
    events <- matrix(c(1, 1, 2, 1), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "event")

    fit <- rem(events = events, formula = ~intercept, n_nodes = 3, directed = TRUE)
    expect_s3_class(fit, "rem")
})

test_that("rem() handles build_time argument through control", {
    events <- matrix(c(
        1, 1, 2, 1,
        2, 2, 3, 1,
        3, 3, 1, 1
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "event")

    fit <- rem(
        events = events,
        formula = ~intercept,
        n_nodes = 4,
        directed = FALSE,
        control = control.redeem(build_time = 1.5)
    )

    expect_s3_class(fit, "rem")
})

test_that("build_time correctly filters and initializes statistics for REM", {
    library(data.table)
    
    # 3 nodes directed sequence:
    # t=1.0: 1 -> 2
    # t=2.0: 2 -> 3
    # t=3.0: 1 -> 2
    # t=3.5: 2 -> 3
    # t=4.0: 3 -> 1
    events <- matrix(c(
        1.0, 1, 2, 1,
        2.0, 2, 3, 1,
        3.0, 1, 2, 1,
        3.5, 2, 3, 1,
        4.0, 3, 1, 1
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "event")

    # Set build_time = 2.5
    fit <- rem(
        events = events,
        formula = ~ intercept + inertia(),
        n_nodes = 3,
        directed = TRUE,
        control = control.redeem(build_time = 2.5, return_data = TRUE)
    )

    # Assess preprocessed data
    prep_data <- as.data.table(fit$model$data)

    # Identifiers
    inertia_col <- grep("inertia", colnames(prep_data), value = TRUE)

    # Check times: Only times >= build_time (2.5) should be present
    expect_true(all(prep_data$time >= 2.5))

    # There should be exactly 2 rows in the preprocessed estimation dataset,
    # both ending at time_new == 4.0:
    # 1. 1 -> 2 starting at time = 3.0, with inertia = 2.
    # 2. 2 -> 3 starting at time = 3.5, with inertia = 2.
    expect_equal(nrow(prep_data), 2)
    expect_true(all(prep_data$time_new == 4.0))

    row_1_2 <- prep_data[from == 1 & to == 2]
    expect_equal(nrow(row_1_2), 1)
    expect_equal(row_1_2$time, 3.0)
    expect_equal(row_1_2[[inertia_col]], 2)

    row_2_3 <- prep_data[from == 2 & to == 3]
    expect_equal(nrow(row_2_3), 1)
    expect_equal(row_2_3$time, 3.5)
    expect_equal(row_2_3[[inertia_col]], 2)
})


