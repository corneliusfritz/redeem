library(testthat)
library(redeem)
library(data.table)

test_that("Participation Shift (P-shift) statistics calculate correctly", {
    # 5 nodes (1, 2, 3, 4, 5)
    # Sequence of directed events:
    # 1. t=1.0: 1 -> 2
    # 2. t=2.0: 2 -> 1 (Reciprocation AB-BA)
    # 3. t=3.0: 1 -> 3 (Turn-Continuing AB-AY)
    # 4. t=4.0: 4 -> 1 (Turn-Usurping AB-XA)
    # 5. t=5.0: 4 -> 2 (Turn-Continuing AB-AY)
    # 6. t=6.0: 3 -> 4 (Full shift AB-XY)
    # 7. t=7.0: 3 -> 4 (Exact repeat/inertia)
    
    events <- matrix(c(
        1.0, 1, 2, 1,
        2.0, 2, 1, 1,
        3.0, 1, 3, 1,
        4.0, 4, 1, 1,
        5.0, 4, 2, 1,
        6.0, 3, 4, 1,
        7.0, 3, 4, 1
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "type")
    
    n_nodes <- 5
    
    # Formula including all 6 P-shifts
    formula <- ~ psABBA() + psABBY() + psABAY() + psABXA() + psABXB() + psABXY()
    
    fit <- rem(
        events = events,
        n_nodes = n_nodes,
        formula = formula,
        directed = TRUE,
        control = control.redeem(return_data = TRUE)
    )
    
    prep_data <- as.data.table(fit$model$data)
    
    # Identify the P-shift columns dynamically to match Rcpp Exports
    psABBA_col <- grep("psABBA", colnames(prep_data), value = TRUE)
    psABBY_col <- grep("psABBY", colnames(prep_data), value = TRUE)
    psABAY_col <- grep("psABAY", colnames(prep_data), value = TRUE)
    psABXA_col <- grep("psABXA", colnames(prep_data), value = TRUE)
    psABXB_col <- grep("psABXB", colnames(prep_data), value = TRUE)
    psABXY_col <- grep("psABXY", colnames(prep_data), value = TRUE)
    
    # 1. Interval [1.0, 2.0): 
    # Previous event was 1 -> 2. So A = 1, B = 2.
    # Expected:
    # - psABBA: active only on 2 -> 1
    # - psABBY: active on 2 -> Y (Y != 1, 2), i.e., 2->3, 2->4, 2->5
    # - psABAY: active on 1 -> Y (Y != 1, 2), i.e., 1->3, 1->4, 1->5
    # - psABXA: active on X -> 1 (X != 1, 2), i.e., 3->1, 4->1, 5->1
    # - psABXB: active on X -> 2 (X != 1, 2), i.e., 3->2, 4->2, 5->2
    # - psABXY: active on X -> Y (X, Y != 1, 2 and X != Y)
    
    int_1_2 <- prep_data[time == 1.0]
    
    # Verify psABBA:
    expect_equal(int_1_2[from == 2 & to == 1, get(psABBA_col)], 1)
    expect_equal(sum(int_1_2[[psABBA_col]]), 1) # Only one dyad is psABBA
    
    # Verify psABBY:
    expect_equal(int_1_2[from == 2 & to %in% c(3, 4, 5), get(psABBY_col)], c(1, 1, 1))
    expect_equal(sum(int_1_2[[psABBY_col]]), 3) # Only 3 dyads
    
    # Verify psABAY:
    expect_equal(int_1_2[from == 1 & to %in% c(3, 4, 5), get(psABAY_col)], c(1, 1, 1))
    expect_equal(sum(int_1_2[[psABAY_col]]), 3)
    
    # Verify psABXA:
    expect_equal(int_1_2[from %in% c(3, 4, 5) & to == 1, get(psABXA_col)], c(1, 1, 1))
    expect_equal(sum(int_1_2[[psABXA_col]]), 3)
    
    # Verify psABXB:
    expect_equal(int_1_2[from %in% c(3, 4, 5) & to == 2, get(psABXB_col)], c(1, 1, 1))
    expect_equal(sum(int_1_2[[psABXB_col]]), 3)
    
    # Verify psABXY:
    expect_equal(int_1_2[!(from %in% c(1, 2)) & !(to %in% c(1, 2)), get(psABXY_col)], rep(1, 6))
    expect_equal(sum(int_1_2[[psABXY_col]]), 6) # 3 * 2 = 6 dyads
    
    # 2. Interval [2.0, 3.0):
    # Previous event was 2 -> 1. So A = 2, B = 1.
    # Expected:
    # - psABBA: active only on 1 -> 2
    # - psABBY: active on 1 -> Y (Y != 1, 2), i.e., 1->3, 1->4, 1->5
    # - psABAY: active on 2 -> Y (Y != 1, 2), i.e., 2->3, 2->4, 2->5
    # - psABXA: active on X -> 2 (X != 1, 2), i.e., 3->2, 4->2, 5->2
    # - psABXB: active on X -> 1 (X != 1, 2), i.e., 3->1, 4->1, 5->1
    # - psABXY: active on X -> Y (X, Y != 1, 2 and X != Y)
    
    int_2_3 <- prep_data[time == 2.0]
    expect_equal(int_2_3[from == 1 & to == 2, get(psABBA_col)], 1)
    expect_equal(sum(int_2_3[[psABBA_col]]), 1)
    expect_equal(int_2_3[from == 1 & to %in% c(3, 4, 5), get(psABBY_col)], c(1, 1, 1))
    expect_equal(int_2_3[from == 2 & to %in% c(3, 4, 5), get(psABAY_col)], c(1, 1, 1))
    expect_equal(int_2_3[from %in% c(3, 4, 5) & to == 2, get(psABXA_col)], c(1, 1, 1))
    expect_equal(int_2_3[from %in% c(3, 4, 5) & to == 1, get(psABXB_col)], c(1, 1, 1))
    expect_equal(sum(int_2_3[[psABXY_col]]), 6)
    
    # 3. Interval [3.0, 4.0):
    # Previous event was 1 -> 3. So A = 1, B = 3.
    # Expected:
    # - psABBA: active only on 3 -> 1
    # - psABBY: active on 3 -> Y (Y != 1, 3), i.e., 3->2, 3->4, 3->5
    # - psABAY: active on 1 -> Y (Y != 1, 3), i.e., 1->2, 1->4, 1->5
    # - psABXA: active on X -> 1 (X != 1, 3), i.e., 2->1, 4->1, 5->1
    # - psABXB: active on X -> 3 (X != 1, 3), i.e., 2->3, 4->3, 5->3
    # - psABXY: active on X -> Y (X, Y != 1, 3 and X != Y)
    
    int_3_4 <- prep_data[time == 3.0]
    expect_equal(int_3_4[from == 3 & to == 1, get(psABBA_col)], 1)
    expect_equal(sum(int_3_4[[psABBA_col]]), 1)
    expect_equal(int_3_4[from == 3 & to %in% c(2, 4, 5), get(psABBY_col)], c(1, 1, 1))
    expect_equal(int_3_4[from == 1 & to %in% c(2, 4, 5), get(psABAY_col)], c(1, 1, 1))
    expect_equal(int_3_4[from %in% c(2, 4, 5) & to == 1, get(psABXA_col)], c(1, 1, 1))
    expect_equal(int_3_4[from %in% c(2, 4, 5) & to == 3, get(psABXB_col)], c(1, 1, 1))
    expect_equal(sum(int_3_4[[psABXY_col]]), 6)
})

test_that("Participation Shift (P-shift) statistics are disallowed for DEM models", {
    events <- matrix(c(
        1.0, 1, 2, 1,
        2.0, 2, 1, 1
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "type")
    
    expect_error(
        dem(
            events = events,
            n_nodes = 5,
            formula_0_1 = ~ psABBA(),
            directed = TRUE
        ),
        "Term 'psABBA' is not allowed for the 'dem' model."
    )
})

test_that("Unified ps(mode) term calculates identical values as specific psABXX terms", {
    events <- matrix(c(
        1.0, 1, 2, 1,
        2.0, 2, 1, 1,
        3.0, 1, 3, 1,
        4.0, 4, 1, 1,
        5.0, 4, 2, 1,
        6.0, 3, 4, 1,
        7.0, 3, 4, 1
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "type")
    
    n_nodes <- 5
    
    # Fit model with individual psABXX terms
    fit_indiv <- rem(
        events = events,
        n_nodes = n_nodes,
        formula = ~ psABBA() + psABBY() + psABAY() + psABXA() + psABXB() + psABXY(),
        directed = TRUE,
        control = control.redeem(return_data = TRUE)
    )
    
    # Fit model with unified ps(mode) terms
    fit_unified <- rem(
        events = events,
        n_nodes = n_nodes,
        formula = ~ ps(mode = "ABBA") + ps(mode = "ABBY") + ps(mode = "ABAY") + ps(mode = "ABXA") + ps(mode = "ABXB") + ps(mode = "ABXY"),
        directed = TRUE,
        control = control.redeem(return_data = TRUE)
    )
    
    # Compare model data matrices
    expect_equal(fit_indiv$model$data, fit_unified$model$data)
    
    # Verify that invalid mode throws an error
    expect_error(
        rem(
            events = events,
            n_nodes = n_nodes,
            formula = ~ ps(mode = "INVALID"),
            directed = TRUE
        ),
        "Mode 'INVALID' for term 'ps' must be one of: ABBA, ABBY, ABAY, ABXA, ABXB, ABXY"
    )
})
