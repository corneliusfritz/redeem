library(testthat)
library(redeem)
library(data.table)

test_that("windowed common partner statistic decays correctly", {
    # Sequence of events to create common partners and then see them expire
    # n_nodes = 4
    # Pairs: (1,3), (2,3), (1,4), (2,4)
    # Window = 5.0
    
    events <- matrix(c(
        1.0, 1, 3, 1,
        2.0, 2, 3, 1,
        3.0, 1, 4, 1,
        4.0, 2, 4, 1,
        10.0, 1, 2, 1 # Dummy event to force calculation at t=10
    ), ncol = 4, byrow = TRUE)
    colnames(events) <- c("time", "from", "to", "type")
    
    n_nodes <- 4
    window_size <- 5
    
    # We use a formula that standardizes to the new naming convention
    formula <- ~ common_partner(window = window_size)
    
    # Fit the model (we don't care about coefficients, just preprocessing)
    fit <- rem(
        events = events,
        n_nodes = n_nodes,
        formula = formula,
        directed = FALSE,
        control = control.redeem(return_data = TRUE)
    )
    
    prep_data <- as.data.table(fit$model$data)
    
    # Identify the common partner column
    stat_col <- grep("common_partner", colnames(prep_data), value = TRUE)
    if (length(stat_col) == 0) {
        cat("\nAvailable columns: ", paste(colnames(prep_data), collapse=", "), "\n")
    }
    expect_length(stat_col, 1)
    stat_col <- stat_col[1]
    
    # Inspect pair (1, 2)
    # Pair ID for (1, 2) in undirected 4-node network:
    # 1-2 is pair 1
    pair1_data <- prep_data[pair_id == 1]
    
    # Let's see the history of pair (1, 2)
    # Intervals are created by all events in the network.
    # Event 1: (1,3) at t=1.0. 
    # Event 2: (2,3) at t=2.0. -> CP(1,2) becomes 1.
    # Event 3: (1,4) at t=3.0.
    # Event 4: (2,4) at t=4.0. -> CP(1,2) becomes 2.
    # Event 5: Expire (1,3) at t=6.0. -> CP(1,2) becomes 1.
    # Event 6: Expire (2,3) at t=7.0. -> CP(1,2) becomes 1? No, 
    #   if 1~3 is gone, then 3 is no longer a common partner.
    #   So at t=6.0, CP(1,2) should drop to 1 (only node 4 is common partner).
    # Event 7: Expire (1,4) at t=8.0. -> CP(1,2) becomes 0.
    
    # Filter for interesting time points
    # We look at the value in each interval [time, time_new)
    
    # --- Manual Calculation for Pair (1, 2) ---
    # t=1.0: (1,3) occurs. Node 3 is a partner of 1.
    # t=2.0: (2,3) occurs. Node 3 is a partner of 2. 
    #        Since 3 is already a partner of 1, (1,2) now has node 3 as a common partner.
    #        CP(1,2) = 1.
    # t=3.0: (1,4) occurs. Node 4 is a partner of 1.
    # t=4.0: (2,4) occurs. Node 4 is a partner of 2.
    #        Since 4 is already a partner of 1, (1,2) now has nodes 3 AND 4 as common partners.
    #        CP(1,2) = 2.
    # t=6.0: Window (size 5) for (1,3) expires (1.0 + 5.0).
    #        Node 3 is no longer a partner of 1.
    #        Therefore, 3 is no longer a common partner for (1,2).
    #        CP(1,2) drops from 2 to 1 (only node 4 remains).
    # t=7.0: Window for (2,3) expires (2.0 + 5.0).
    #        Node 3 is no longer a partner of 2.
    #        (3 was already not a partner of 1, so CP(1,2) remains 1).
    # t=8.0: Window for (1,4) expires (3.0 + 5.0).
    #        Node 4 is no longer a partner of 1.
    #        CP(1,2) drops from 1 to 0.
    
    # Interval [2.0, 3.0): CP should be 1
    expect_equal(pair1_data[time == 2.0, get(stat_col)], 1)
    
    # Interval [4.0, 6.0): CP should be 2
    expect_equal(pair1_data[time == 4.0, get(stat_col)], 2)
    
    # Interval [6.0, 7.0): CP should be 1
    expect_equal(pair1_data[time == 6.0, get(stat_col)], 1)
    
    # Interval [8.0, 9.0): CP should be 0
    expect_equal(pair1_data[time == 8.0, get(stat_col)], 0)
    
    # Check that it actually went down
    diffs <- diff(pair1_data[[stat_col]])
    expect_true(any(diffs < 0))
})
