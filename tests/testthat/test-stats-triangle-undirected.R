library(testthat)
library(redeem)

test_that("Undirected triangle statistics calculation", {
    n_nodes <- 3
    # Matrix must have 4 columns for preprocess_rem: time, from, to, type
    events <- matrix(c(
        1.0, 1, 3, 1,
        2.0, 2, 3, 1,
        3.0, 1, 2, 1,
        4.0, 1, 2, 1
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "type")

    pre <- redeem:::formula_preprocess(
        formula_0_1 = ~ general_triangle(),
        formula_1_0 = ~ general_triangle(),
        events = events,
        n_nodes = n_nodes,
        model_type = "rem"
    )

    out <- redeem:::preprocess_rem(
        edgelist = as.matrix(pre$events),
        n_nodes = n_nodes,
        terms = pre$term_names,
        data_list = pre$data_list,
        transformations = pre$transformation_list,
        directed = FALSE,
        verbose = FALSE
    )

    df <- as.data.frame(out)
    colnames(df) <- c("time_end", "time_start", "pair_id", "status", "event", "from", "to", "f_av", "t_av", "Triangle")

    # The interval starting at t=3.0 (when 1-2 forms)
    row3 <- df[df$time_start == 3.0 & df$from == 1 & df$to == 2, ]
    expect_equal(nrow(row3), 1)
    expect_equal(row3$Triangle, 1.0)
})

test_that("Undirected triangle updates other dyads", {
    n_nodes <- 3
    events <- matrix(c(
        1.0, 1, 2, 1,
        2.0, 1, 3, 1,
        3.0, 2, 3, 1,
        4.0, 2, 3, 1 # Dummy event to see state after t=3.0
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "type")

    pre <- redeem:::formula_preprocess(
        formula_0_1 = ~ general_triangle(),
        events = events,
        n_nodes = n_nodes,
        model_type = "rem"
    )

    out <- redeem:::preprocess_rem(
        edgelist = as.matrix(pre$events),
        n_nodes = n_nodes,
        terms = pre$term_names,
        data_list = pre$data_list,
        transformations = pre$transformation_list,
        directed = FALSE,
        verbose = FALSE
    )

    df <- as.data.frame(out)
    colnames(df) <- c("time_end", "time_start", "pair_id", "status", "event", "from", "to", "f_av", "t_av", "Triangle")
    # print(df)

    row3_12 <- df[df$time_start == 3.0 & df$from == 1 & df$to == 2, ]
    expect_equal(nrow(row3_12), 1)
    expect_equal(row3_12$Triangle, 1.0)

    # And focal dyad (2,3) itself at t=3.0
    row3_23 <- df[df$time_start == 3.0 & df$from == 2 & df$to == 3, ]
    expect_equal(nrow(row3_23), 1)
    expect_equal(row3_23$Triangle, 1.0)
})
