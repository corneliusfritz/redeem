library(testthat)
library(redeem)

test_that("degree dispatcher returns identical results to explicit terms", {
    n_nodes <- 3
    events <- data.frame(
        time = c(1.0, 2.0, 3.0),
        from = c(1, 1, 1),
        to = c(2, 3, 2),
        type = c(TRUE, TRUE, FALSE)
    )

    # 1. Test current_degree_sum (Undirected)
    pre_explicit <- redeem:::formula_preprocess(
        formula_0_1 = ~ current_degree_sum(),
        events = events,
        n_nodes = n_nodes,
        model_type = "dem",
        directed = FALSE
    )

    pre_dispatched <- redeem:::formula_preprocess(
        formula_0_1 = ~ degree(history = "current", type = "sum", count = FALSE),
        events = events,
        n_nodes = n_nodes,
        model_type = "dem",
        directed = FALSE
    )

    expect_equal(unname(pre_dispatched$term_names), unname(pre_explicit$term_names))
    expect_equal(names(pre_dispatched$data_list), names(pre_explicit$data_list))
})

test_that("count dispatcher returns identical results to explicit terms", {
    n_nodes <- 3
    events <- data.frame(
        time = c(1.0, 2.0, 3.0),
        from = c(1, 1, 1),
        to = c(2, 3, 2),
        type = c(TRUE, TRUE, FALSE)
    )

    pre_explicit <- redeem:::formula_preprocess(
        formula_0_1 = ~ general_count_out_sender(),
        events = events,
        n_nodes = n_nodes,
        model_type = "rem",
        directed = TRUE
    )

    pre_dispatched <- redeem:::formula_preprocess(
        formula_0_1 = ~ count(history = "general", type = "out_sender"),
        events = events,
        n_nodes = n_nodes,
        model_type = "rem",
        directed = TRUE
    )

    expect_equal(unname(pre_dispatched$term_names), unname(pre_explicit$term_names))
    expect_equal(names(pre_dispatched$data_list), names(pre_explicit$data_list))
    expect_true("general_count_out_sender" %in% pre_dispatched$term_names)
})

test_that("triangle and common_partner dispatchers work", {
    n_nodes <- 3
    events <- data.frame(
        time = c(1.0, 2.0, 3.0),
        from = c(1, 1, 1),
        to = c(2, 3, 2),
        type = c(TRUE, TRUE, FALSE)
    )

    # Triangle (Default type is OSP)
    pre_t_dispatched <- redeem:::formula_preprocess(
        formula_0_1 = ~ triangle(history = "general"),
        events = events,
        n_nodes = n_nodes,
        model_type = "dem",
        directed = TRUE
    )
    expect_true("general_triangle_OSP" %in% pre_t_dispatched$term_names)

    # Common Partner
    pre_cp_dispatched <- redeem:::formula_preprocess(
        formula_0_1 = ~ common_partner(history = "current"),
        events = events,
        n_nodes = n_nodes,
        model_type = "dem",
        directed = TRUE
    )
    expect_true("current_common_partner_OSP" %in% pre_cp_dispatched$term_names)

    # Triangle with specific type (Directed)
    pre_t_isp <- redeem:::formula_preprocess(
        formula_0_1 = ~ triangle(history = "general", type = "ISP"),
        events = events,
        n_nodes = n_nodes,
        model_type = "rem",
        directed = TRUE
    )
    expect_true("general_triangle_ISP" %in% pre_t_isp$term_names)
})

test_that("transformations still add suffixes but not for identity", {
    n_nodes <- 3
    events <- data.frame(time=1, from=1, to=2, type=TRUE)
    
    pre <- redeem:::formula_preprocess(
        formula_0_1 = ~ degree(history="current", type="sum", transformation="log"),
        events = events,
        n_nodes = n_nodes,
        model_type = "dem",
        directed = FALSE
    )
    expect_true("current_degree_sum_log" %in% names(pre$data_list))
    expect_false(any(grepl("identity", names(pre$data_list))))
})
