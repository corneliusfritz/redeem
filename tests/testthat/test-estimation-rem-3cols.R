library(testthat)
library(redeem)

test_that("rem() works with 3-column event matrix", {
    n <- 5
    # Only time, from, to
    events <- matrix(c(
        1.0, 1, 3,
        2.0, 3, 2,
        3.0, 1, 2
    ), ncol = 3, byrow = TRUE)
    
    # This should not error
    expect_error(
        fit <- rem(
            events = events,
            n_nodes = n,
            formula = ~ 1,
            control = control.redeem(it_max = 5)
        ),
        NA
    )
    
    expect_s3_class(fit, "rem")
    expect_equal(ncol(fit$events), 4)
    expect_true(all(fit$events[, 4] == 1))
})

test_that("dem() provides helpful error for 3-column matrix", {
    n <- 5
    events <- matrix(c(
        1.0, 1, 3,
        2.0, 3, 2
    ), ncol = 3, byrow = TRUE)
    
    expect_error(
        dem(
            events = events,
            n_nodes = n,
            formula_0_1 = ~ 1,
            formula_1_0 = ~ 1
        ),
        "Event matrix must have at least 4 columns"
    )
})
