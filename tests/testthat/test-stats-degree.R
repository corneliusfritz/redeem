library(testthat)
library(redeem)

test_that("Degree statistics (General/Current, Degree/Count) calculation", {
    # Small hand-checkable example
    # Nodes 1, 2, 3
    # t=1: 1 -> 2 (Formation)
    # t=2: 1 -> 3 (Formation)
    # t=3: 1 -> 2 (Dissolution)
    # t=4: 1 -> 2 (Formation)
    
    n_nodes <- 3
    events <- matrix(c(
        1.0, 1, 2, 1,
        2.0, 1, 3, 1,
        3.0, 1, 2, 0,
        4.0, 1, 2, 1
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "type")

    # row at t=0.0 (Interval [0, 1]): State before first event
    #   gen_deg=0, gen_cnt=0, cur_deg=0, cur_cnt=0.
    # row at t=1.0 (Interval [1, 2]): State after 1->2 forms. 
    #   gen_deg=1 (2), gen_cnt=1, cur_deg=1, cur_cnt=1.
    # row at t=2.0 (Interval [2, 3]): State after 1->3 forms.
    #   gen_deg=2 (2,3), gen_cnt=2, cur_deg=2, cur_cnt=2.
    # row at t=3.0 (Interval [3, 4]): State after 1->2 dissolves.
    #   gen_deg=2, gen_cnt=2, cur_deg=1 (only 3), cur_cnt=1.

    pre <- redeem:::formula_preprocess(
        formula_0_1 = ~ general_degree_out_sender() + general_count_out_sender() + 
                        current_degree_out_sender() + current_count_out_sender(),
        formula_1_0 = ~ general_degree_out_sender() + general_count_out_sender() + 
                        current_degree_out_sender() + current_count_out_sender(),
        events = events,
        n_nodes = n_nodes,
        model_type = "dem",
        directed = TRUE
    )

    out <- redeem:::preprocess(
        edgelist = as.matrix(pre$events),
        terms = pre$term_names,
        data_list = pre$data_list,
        transformations = pre$transformation_list,
        n_nodes = n_nodes,
        verbose = FALSE,
        directed = TRUE,
        simultaneous_interactions = FALSE
    )

    df <- as.data.frame(out)
    colnames(df) <- c("time_end", "time_start", "pair_id", "status", "event", "from", "to", "f_av", "t_av", 
                      "gen_deg", "gen_cnt", "cur_deg", "cur_cnt")

    # Before any event (Interval [0, 1])
    row0 <- df[df$time_start == 0 & df$from == 1 & df$to == 2, ]
    expect_equal(row0$gen_deg, 0)
    expect_equal(row0$cur_deg, 0)

    # After t=1.0 formation (Interval [1, 2])
    row1 <- df[df$time_start == 1.0 & df$from == 1 & df$to == 2, ]
    expect_equal(row1$gen_deg, 1)
    expect_equal(row1$gen_cnt, 1)
    expect_equal(row1$cur_deg, 1)
    expect_equal(row1$cur_cnt, 1)

    # After t=2.0 formation (Interval [2, 3])
    row2 <- df[df$time_start == 2.0 & df$from == 1 & df$to == 3, ]
    expect_equal(row2$gen_deg, 2)
    expect_equal(row2$gen_cnt, 2)
    expect_equal(row2$cur_deg, 2)
    expect_equal(row2$cur_cnt, 2)

    # After t=3.0 dissolution (Interval [3, 4])
    row3 <- df[df$time_start == 3.0 & df$from == 1 & df$to == 2, ]
    expect_equal(row3$gen_deg, 2)
    expect_equal(row3$gen_cnt, 2)
    expect_equal(row3$cur_deg, 1)
    expect_equal(row3$cur_cnt, 1)
})

test_that("Undirected Degree statistics calculation", {
    n_nodes <- 3
    # Add events: 1-2 forms (t=1), 1-3 forms (t=2), 1-2 dissolves (t=3)
    # Add a dummy event at 3.1 to see the state AFTER t=3.0
    events <- data.frame(
        time = c(1.0, 2.0, 3.0, 3.1),
        from = c(1, 1, 1, 2),
        to = c(2, 3, 2, 3),
        type = c(TRUE, TRUE, FALSE, TRUE)
    )

    pre <- redeem:::formula_preprocess(
        formula_0_1 = ~ general_degree_sum() + general_count_sum() + 
                        current_degree_sum() + current_count_sum(),
        events = events,
        n_nodes = n_nodes,
        model_type = "dem",
        directed = FALSE
    )

    out <- redeem:::preprocess(
        edgelist = as.matrix(pre$events),
        terms = pre$term_names,
        data_list = pre$data_list,
        transformations = pre$transformation_list,
        n_nodes = n_nodes,
        verbose = FALSE,
        directed = FALSE,
        simultaneous_interactions = FALSE
    )

    df <- as.data.frame(out)
    colnames(df) <- c("time_end", "time_start", "pair_id", "status", "event", "from", "to", "f_av", "t_av", 
                      "gen_deg_sum", "gen_cnt_sum", "cur_deg_sum", "cur_cnt_sum")
    
    # t=1.0: 1-2 forms. Sum=2.
    row1 <- df[abs(df$time_start - 1.0) < 1e-7 & df$from == 1 & df$to == 2, ]
    expect_equal(row1$gen_deg_sum, 2)

    # t=2.0: 1-3 forms. Sum=3.
    row2 <- df[abs(df$time_start - 2.0) < 1e-7 & df$from == 1 & df$to == 3, ]
    expect_equal(row2$gen_deg_sum, 3)
    expect_equal(row2$cur_deg_sum, 3)

    # t=3.0: 1-2 dissolves. Cur Sum=1 (only 1-3 remains).
    row3 <- df[abs(df$time_start - 3.0) < 1e-7, ]
    row3_23 <- row3[row3$from == 2 & row3$to == 3, ]
    expect_equal(row3_23$cur_deg_sum, 1)
})
