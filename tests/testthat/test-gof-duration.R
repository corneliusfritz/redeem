library(testthat)
library(redeem)

test_that("get_ranking() works correctly with current_interaction", {
    set.seed(123)
    n_nodes <- 5
    
    # Simulate data with current_interaction (duration effect)
    events <- dem.simulate(
        formula_0_1 = ~ intercept,
        formula_1_0 = ~ intercept + current_interaction,
        coef_0_1 = c("intercept" = -3),
        coef_1_0 = c("intercept" = -2, "current_interaction" = 0.5),
        n_events = 50,
        n_nodes = n_nodes,
        directed = TRUE
    )
    
    fit <- dem(
        events = events[1:40, ],
        formula_0_1 = ~ intercept,
        formula_1_0 = ~ intercept + current_interaction,
        n_nodes = n_nodes,
        directed = TRUE
    )
    
    # Ranking on the last 10 events
    ranking <- get_ranking(fit, edgelist_test = events[41:50, ], k_max = 10)
    
    expect_s3_class(ranking, "data.frame")
    expect_true(all(ranking$Recall >= 0))
})
