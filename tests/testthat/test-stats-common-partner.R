library(testthat)
library(redeem)

test_that("directed common partner dispatchers work", {
    n_nodes <- 3
    events <- data.frame(
        time = c(1.0, 2.0),
        from = c(1, 2),
        to = c(2, 3),
        type = c(TRUE, TRUE)
    )

    # ISP
    pre_isp <- redeem:::formula_preprocess(
        formula_0_1 = ~ common_partner(history="general", type="ISP"),
        events = events,
        n_nodes = n_nodes,
        model_type = "rem",
        directed = TRUE
    )
    expect_true("general_common_partner_ISP" %in% pre_isp$term_names)

    # OSP
    pre_osp <- redeem:::formula_preprocess(
        formula_0_1 = ~ common_partner(history="general", type="OSP"),
        events = events,
        n_nodes = n_nodes,
        model_type = "rem",
        directed = TRUE
    )
    expect_true("general_common_partner_OSP" %in% pre_osp$term_names)

    # OTP
    pre_otp <- redeem:::formula_preprocess(
        formula_0_1 = ~ common_partner(history="general", type="OTP"),
        events = events,
        n_nodes = n_nodes,
        model_type = "rem",
        directed = TRUE
    )
    expect_true("general_common_partner_OTP" %in% pre_otp$term_names)

    # ITP
    pre_itp <- redeem:::formula_preprocess(
        formula_0_1 = ~ common_partner(history="general", type="ITP"),
        events = events,
        n_nodes = n_nodes,
        model_type = "rem",
        directed = TRUE
    )
    expect_true("general_common_partner_ITP" %in% pre_itp$term_names)
})

test_that("undirected common partner still works", {
    n_nodes <- 3
    events <- data.frame(time=1, from=1, to=2, type=TRUE)
    
    pre <- redeem:::formula_preprocess(
        formula_0_1 = ~ common_partner(history="general"),
        events = events,
        n_nodes = n_nodes,
        model_type = "dem",
        directed = FALSE
    )
    expect_true("general_common_partner" %in% pre$term_names)
    expect_false(any(grepl("OSP", pre$term_names)))
})
