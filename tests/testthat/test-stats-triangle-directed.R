library(testthat)
library(redeem)

get_stat_val <- function(df, from_node, to_node, t, stat_col) {
    sub <- df[df$from == from_node & df$to == to_node & df$time_start <= t, ]
    if (nrow(sub) == 0) return(0.0)
    sub[which.max(sub$time_start), stat_col]
}

test_that("Directed triangle statistics (OSP) calculation", {
    # OSP(i, j): i -> k, j -> k
    # Closed triangle: statistic is 1.0 ONLY if i -> j also exists.
    n_nodes <- 3
    events <- matrix(c(
        1.0, 2, 3, 1,
        2.0, 1, 3, 1,
        3.0, 1, 2, 1,
        4.0, 1, 2, 1 # Dummy
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "type")

    pre <- redeem:::formula_preprocess(
        formula_0_1 = ~ general_triangle(type = "OSP"),
        formula_1_0 = ~ general_triangle(type = "OSP"),
        events = events,
        n_nodes = n_nodes,
        model_type = "rem", directed = TRUE
    )

    out <- redeem:::preprocess_rem(
        edgelist = as.matrix(pre$events),
        n_nodes = n_nodes,
        terms = pre$term_names,
        data_list = pre$data_list,
        transformations = pre$transformation_list,
        directed = TRUE,
        verbose = FALSE
    )

    df <- as.data.frame(out)
    colnames(df) <- c("time_end", "time_start", "pair_id", "status", "event", "from", "to", "f_av", "t_av", "OSP")
    row3 <- df[df$time_start == 3.0 & df$from == 1 & df$to == 2, ]
    expect_equal(nrow(row3), 1)
    expect_equal(row3$OSP, 1.0)

    # OSP(2,1) is still 0.0 because 2->1 has NOT formed
    row3_rev <- df[df$time_start == 0.0 & df$from == 2 & df$to == 1, ]
    expect_equal(nrow(row3_rev), 1)
    expect_equal(row3_rev$OSP, 0.0)
})

test_that("Directed triangle statistics (OTP) calculation", {
    # OTP(i, j): i -> k -> j
    n_nodes <- 3
    events <- matrix(c(
        1.0, 2, 3, 1,
        2.0, 1, 2, 1,
        3.0, 1, 3, 1,
        4.0, 1, 3, 1 # Dummy
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "type")

    pre <- redeem:::formula_preprocess(
        formula_0_1 = ~ general_triangle(type = "OTP"),
        formula_1_0 = ~ general_triangle(type = "OTP"),
        events = events,
        n_nodes = n_nodes,
        model_type = "rem", directed = TRUE
    )

    out <- redeem:::preprocess_rem(
        edgelist = as.matrix(pre$events),
        n_nodes = n_nodes,
        terms = pre$term_names,
        data_list = pre$data_list,
        transformations = pre$transformation_list,
        directed = TRUE,
        verbose = FALSE
    )
    df <- as.data.frame(out)
    colnames(df) <- c("time_end", "time_start", "pair_id", "status", "event", "from", "to", "f_av", "t_av", "OTP")

    # At t=3.0, (1, 3) formed. 1->2 exists, 2->3 exists. So OTP(1,3) becomes 1.0
    row3 <- df[df$time_start == 3.0 & df$from == 1 & df$to == 3, ]
    expect_equal(nrow(row3), 1)
    expect_equal(row3$OTP, 1.0)
})


test_that("Asymmetry of OTP directed triangle statistics", {
    # 3 nodes:
    # t=1.0: 1 -> 2
    # t=2.0: 2 -> 3
    # t=3.0: 1 -> 3 (transitive OTP: 1 -> 2 -> 3)
    # t=4.0: 1 -> 2 (dummy to keep t=3.0 interval)
    n_nodes <- 3
    events <- matrix(c(
        1.0, 1, 2, 1,
        2.0, 2, 3, 1,
        3.0, 1, 3, 1,
        4.0, 1, 2, 1
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "type")

    pre <- redeem:::formula_preprocess(
        formula_0_1 = ~ general_triangle(type = "OTP"),
        events = events,
        n_nodes = n_nodes,
        model_type = "rem", directed = TRUE
    )

    out <- redeem:::preprocess_rem(
        edgelist = as.matrix(pre$events),
        n_nodes = n_nodes,
        terms = pre$term_names,
        data_list = pre$data_list,
        transformations = pre$transformation_list,
        directed = TRUE,
        verbose = FALSE
    )
    df <- as.data.frame(out)
    colnames(df) <- c("time_end", "time_start", "pair_id", "status", "event", "from", "to", "f_av", "t_av", "OTP")

    # Before t=3.0, (1,3) has not occurred/been observed, so OTP = 0.0
    expect_equal(get_stat_val(df, 1, 3, 2.0, "OTP"), 0.0)

    # After t=3.0, (1,3) has occurred/been observed, so OTP = 1.0
    expect_equal(get_stat_val(df, 1, 3, 3.0, "OTP"), 1.0)

    # Reverse direction (3 -> 1) should have OTP = 0.0 (asymmetry check)
    expect_equal(get_stat_val(df, 3, 1, 3.0, "OTP"), 0.0)
})

test_that("Asymmetry of ITP directed triangle statistics", {
    # 3 nodes:
    # t=1.0: 3 -> 2
    # t=2.0: 2 -> 1
    # t=3.0: 1 -> 3 (transitive ITP: 3 -> 2 -> 1)
    # t=4.0: 1 -> 2 (dummy to keep t=3.0 interval)
    n_nodes <- 3
    events <- matrix(c(
        1.0, 3, 2, 1,
        2.0, 2, 1, 1,
        3.0, 1, 3, 1,
        4.0, 1, 2, 1
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "type")

    pre <- redeem:::formula_preprocess(
        formula_0_1 = ~ general_triangle(type = "ITP"),
        events = events,
        n_nodes = n_nodes,
        model_type = "rem", directed = TRUE
    )

    out <- redeem:::preprocess_rem(
        edgelist = as.matrix(pre$events),
        n_nodes = n_nodes,
        terms = pre$term_names,
        data_list = pre$data_list,
        transformations = pre$transformation_list,
        directed = TRUE,
        verbose = FALSE
    )
    df <- as.data.frame(out)
    colnames(df) <- c("time_end", "time_start", "pair_id", "status", "event", "from", "to", "f_av", "t_av", "ITP")

    # Before t=3.0, (1,3) has not occurred/been observed, so ITP = 0.0
    expect_equal(get_stat_val(df, 1, 3, 2.0, "ITP"), 0.0)

    # After t=3.0, (1,3) has occurred/been observed, so ITP = 1.0 (since 3 -> 2 -> 1 exists)
    expect_equal(get_stat_val(df, 1, 3, 3.0, "ITP"), 1.0)

    # Reverse direction (3 -> 1) should have ITP = 0.0 (asymmetry check)
    expect_equal(get_stat_val(df, 3, 1, 3.0, "ITP"), 0.0)
})

test_that("Comparison of common_partners and triangle statistics", {
    # 3 nodes:
    # t=1.0: 1 -> 2
    # t=2.0: 3 -> 2
    # t=3.0: 1 -> 3
    # t=4.0: 1 -> 3 (dummy)
    n_nodes <- 3
    events <- matrix(c(
        1.0, 1, 2, 1,
        2.0, 3, 2, 1,
        3.0, 1, 3, 1,
        4.0, 1, 3, 1
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "type")

    # In this network:
    # At t=2.0 (interval [2.0, 3.0]), 1->2 and 3->2 exist, sharing receiver 2.
    # Therefore, OSP shared partner for (1,3) is 2.
    # Since 1->3 has not yet formed, common_partners(type="OSP") for (1,3) should be 1.0,
    # but triangle(type="OSP") for (1,3) should be 0.0 (since the edge 1->3 is inactive/unobserved).
    # At t=3.0, 1->3 is formed. Now triangle(type="OSP") for (1,3) should become 1.0.

    pre <- redeem:::formula_preprocess(
        formula_0_1 = ~ general_common_partners(type = "OSP") + general_triangle(type = "OSP"),
        events = events,
        n_nodes = n_nodes,
        model_type = "rem", directed = TRUE
    )

    out <- redeem:::preprocess_rem(
        edgelist = as.matrix(pre$events),
        n_nodes = n_nodes,
        terms = pre$term_names,
        data_list = pre$data_list,
        transformations = pre$transformation_list,
        directed = TRUE,
        verbose = FALSE
    )
    df <- as.data.frame(out)
    colnames(df) <- c("time_end", "time_start", "pair_id", "status", "event", "from", "to", "f_av", "t_av", "CP_OSP", "Tri_OSP")

    # At t=2.0, CP_OSP(1,3) is 1.0, Tri_OSP(1,3) is 0.0
    expect_equal(get_stat_val(df, 1, 3, 2.0, "CP_OSP"), 1.0)
    expect_equal(get_stat_val(df, 1, 3, 2.0, "Tri_OSP"), 0.0)

    # At t=3.0, CP_OSP(1,3) is 1.0, Tri_OSP(1,3) is 1.0
    expect_equal(get_stat_val(df, 1, 3, 3.0, "CP_OSP"), 1.0)
    expect_equal(get_stat_val(df, 1, 3, 3.0, "Tri_OSP"), 1.0)
})

