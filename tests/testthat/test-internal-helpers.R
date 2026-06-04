library(testthat)
library(redeem)

test_that("formula_preprocess_single() works", {
    # Basic formula
    f <- ~ 1
    p <- redeem:::formula_preprocess_single(f, n_nodes = 5)
    expect_equal(unname(p$coef_names), "Intercept")
    expect_false(p$includes_degrees)

    # Formula with degrees
    f_deg <- ~ degrees + dyad_cov
    p_deg <- redeem:::formula_preprocess_single(f_deg, n_nodes = 5)
    expect_true(p_deg$includes_degrees)
})

test_that("formula_preprocess() combines correctly", {
    events <- matrix(c(0,0,0,0), nrow=1)
    colnames(events) <- c("time", "from", "to", "event")

    p <- redeem:::formula_preprocess(
        formula_0_1 = ~ 1,
        formula_1_0 = ~ 1,
        events = events,
        n_nodes = 5
    )

    expect_true("Intercept" %in% p$coef_names)
    expect_equal(p$n_nodes, 5)
})

test_that("rhs_terms_as_list() evaluates symbols", {
    cov_val <- matrix(1, 3, 3)
    f <- ~ cov_val
    # Note: rhs_terms_as_list is internal, but we can test it via preprocess_single or
    # if it's exported for testing. It's not exported, so we test it implicitly.
    p <- redeem:::formula_preprocess_single(f, n_nodes = 3)
    expect_true("cov_val" %in% p$coef_names)
})
