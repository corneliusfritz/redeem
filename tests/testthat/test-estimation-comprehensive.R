
library(testthat)
library(redeem)
library(data.table)

test_that("Directed Statistics (Degrees, Triadic, Reciprocity) are correct by hand", {
  n_nodes <- 3
  
  # Event Sequence:
  # t=1: 1 -> 2 (F)
  # t=2: 2 -> 3 (F)
  # t=3: 1 -> 3 (F) -> Trigger point for triangles
  # t=4: 1 -> 2 (D) -> Breaking triadic chain
  # t=5: 1 -> 3 (D) -> Breaking last link
  
  events <- matrix(c(
    1.0, 1, 2, 1,
    2.0, 2, 3, 1,
    3.0, 1, 3, 1,
    4.0, 1, 2, 0,
    5.0, 1, 3, 0,
    10.0, 2, 3, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")
  
  # Statistics to calculate:
  # 1. inertia()
  # 2. current_triangle(type = "OTP")
  # 3. general_triangle(type = "OTP")
  # 4. current_degree_out_sender()
  # 5. general_degree_out_sender()
  # 6. reciprocity()
  
  formula_0_1 <- ~ Intercept + 
                   inertia() + 
                   current_triangle(type = "OTP") + 
                   general_triangle(type = "OTP") + 
                   current_degree_out_sender() + 
                   general_degree_out_sender() + 
                   reciprocity()
                   
  fit <- dem(
    events = events,
    formula_0_1 = formula_0_1,
    formula_1_0 = ~ Intercept,
    n_nodes = n_nodes,
    directed = TRUE,
    control = control.redeem(return_data = TRUE)
  )
  
  # Asses returned data
  # "Intercept", "inertia_identity", "current_triangle_OTP_identity", etc.
  dt <- as.data.table(fit$model_0_1$data)
  
  # Dyad (1, 2) dissolved at t=4, so it should be in model_0_1 after t=4.
  dyad_1_2 <- dt[from == 1 & to == 2]
  dyad_1_2 <- dyad_1_2[order(time)]
  
  # Interval [4, 5]: After 1->2 dissolution.
  # OTP_curr=0 (chain broken), OTP_gen=0 (since 1->3 formed while 1->2 was active? no, 1->2 focal).
  # Wait, OTP for (1,2) involves (1->k, k->2). At t=4, 1->3 is active. But 3->2? No.
  # general_degree_out(1) should be 1 (because 1->3 is STILL active? No, general_ is history).
  # 1->2 had 1 event. 1->3 had 1 event. So gen_deg_out(1) = 2.
  row_t4_5 <- dyad_1_2[time == 4]
  expect_equal(row_t4_5$current_triangle_OTP_identity, 0)
  expect_equal(row_t4_5$current_degree_out_sender_identity, 1) # 1->3 is currently active
  expect_equal(row_t4_5$general_degree_out_sender_identity, 2) # Both 1->2 and 1->3 happened
  
  message("Directed Statistics Correct.")
})


test_that("Undirected Statistics (Common Partner, Degree Sum) are correct by hand", {
  n_nodes <- 3
  
  # Event Sequence:
  # t=1: 1 - 2 (F)
  # t=2: 2 - 3 (F)
  # t=4: 1 - 2 (D)
  # t=5: 2 - 3 (D) [Optional cleanup]
  
  events <- matrix(c(
    1.0, 1, 2, 1,
    2.0, 2, 3, 1,
    4.0, 1, 2, 0,
    5.0, 2, 3, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")
  
  formula_0_1 <- ~ Intercept + 
                   current_common_partners() + 
                   general_common_partners() + 
                   degree_sum() + 
                   degree_absdiff()
                   
  fit <- dem(
    events = events,
    formula_0_1 = formula_0_1,
    formula_1_0 = ~ Intercept,
    n_nodes = n_nodes,
    directed = FALSE,
    control = control.redeem(return_data = TRUE)
  )
  
  dt <- as.data.table(fit$model_0_1$data)
  
  # Dyad (1, 3) in undirected mode (from < to)
  dyad_1_3 <- dt[from == 1 & to == 3]
  dyad_1_3 <- dyad_1_3[order(time)]
  
  col_curr <- grep("current_common_partner", names(dt), value = TRUE)
  col_gen <- grep("general_common_partner", names(dt), value = TRUE)
  col_sum <- grep("degree_sum", names(dt), value = TRUE)
  
  # Interval [2, 4]: After 1-2 and 2-3 formed.
  # (1, 3) common partner: 1 (node 2).
  # degree_sum: deg(1)=1, deg(3)=1 -> sum=2.
  row_t2_4 <- dyad_1_3[time == 2]
  expect_equal(row_t2_4[[col_curr[1]]], 1)
  expect_equal(row_t2_4[[col_gen[1]]], 1)
  expect_equal(row_t2_4[[col_sum[1]]], 2)
  
  # Interval [4, 5]: After 1-2 dissolves.
  # (1, 3) current CP: 0.
  # (1, 3) general CP: 1 (history persists).
  row_t4_5 <- dyad_1_3[time == 4]
  expect_equal(row_t4_5[[col_curr[1]]], 0)
  expect_equal(row_t4_5[[col_gen[1]]], 1)
  
  message("Undirected Statistics Correct.")
})


test_that("Multi-Stream Covariate (Inertia from exogenous stream) is correct", {
  n_nodes <- 3
  
  # Main events (empty timeline mostly)
  # We just need some intervals to check.
  events_main <- matrix(c(
    10.0, 1, 2, 1,
    20.0, 1, 2, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events_main) <- c("time", "from", "to", "type")
  
  # Exogenous stream
  events_exo <- matrix(c(
    1.0, 1, 3, 1,
    5.0, 1, 3, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events_exo) <- c("time", "from", "to", "type")
  
  formula_0_1 <- ~ Intercept + inertia(event_stream = events_exo)
  
  fit <- dem(
    events = events_main,
    formula_0_1 = formula_0_1,
    formula_1_0 = ~ Intercept,
    n_nodes = n_nodes,
    directed = TRUE,
    simultaneous_interactions = TRUE,
    control = control.redeem(
      return_data = TRUE
    )
  )
  
  dt <- as.data.table(fit$model_0_1$data)
  
  # Dyad (1, 3) at main event intervals.
  # Main event at t=10. Risk set intervals include [0, 10].
  # Wait, dem splits [0, 10] into sub-intervals based on exo events?
  # Yes, preprocess_multi_stream splinters the timeline.
  # Sub-intervals for (1,3) should be:
  # [0, 1]: No exo event -> stat 0
  # [1, 5]: Exo event active -> stat 1
  # [5, 10]: Exo event ended -> stat 1 (historical inertia)
  
  dyad_1_3 <- dt[from == 1 & to == 3]
  dyad_1_3 <- dyad_1_3[order(time)]
  print("Multi-Stream Dyad (1,3) intervals:")
  print(dyad_1_3)
  
  # Column name will be inertia.event_stream... 
  stat_col <- grep("inertia", names(dt), value = TRUE)
  print(paste("Stat column found:", stat_col))
  
  expect_equal(dyad_1_3[time == 1, get(stat_col[1])], 0)
  expect_equal(dyad_1_3[time == 5, get(stat_col[1])], 1)
  expect_equal(dyad_1_3[time == 0, get(stat_col[1])], 0)
  
  message("Multi-Stream Statistics Correct.")
})


test_that("Structural Statistics (Intercept, current_interaction) are correct", {
  n_nodes <- 3
  
  # Event Sequence for 1->0 (Dissolution)
  # t=1: 1 -> 2 (F)
  # t=5: 1 -> 2 (D)
  # t=10: dummy end
  
  events <- matrix(c(
    1.0, 1, 2, 1,
    5.0, 1, 2, 0,
    10.0, 1, 2, 1, # Dummy to have an interval after
    20.0, 1, 2, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")
  
  # For 1->0, we can use current_interaction()
  # It measures duration since last event.
  formula_1_0 <- ~ Intercept + current_interaction()
  
  fit <- dem(
    events = events,
    formula_0_1 = ~ Intercept,
    formula_1_0 = formula_1_0,
    n_nodes = n_nodes,
    directed = TRUE,
    control = control.redeem(return_data = TRUE)
  )
  
  dt <- as.data.table(fit$model_1_0$data)
  
  # Dyad (1, 2) is active in [1, 5].
  dyad_1_2 <- dt[from == 1 & to == 2]
  dyad_1_2 <- dyad_1_2[order(time)]
  
  # In interval [1, 5], 1->2 is at risk of dissolution.
  # current_interaction starts at 0 at t=1 and increases?
  # Wait, stat_current_interaction adds (time_now - time_last).
  # During preprocessing, for a long interval it might be calculated at checkpoints.
  # Let's check the value at the beginning of the interval.
  
  # At t=1, it should be 0.
  expect_equal(dyad_1_2[time == 1, current_interaction_identity], 0)
  
  # Intercept should always be 1
  expect_true(all(dt$Intercept == 1))
  
  message("Structural Statistics Correct.")
})

test_that("build_time argument correctly separates burn-in history from the estimation phase", {
  n_nodes <- 3
  
  # Expanded sequence:
  # Dyad (1,2) forms at 1.0 (inside burn-in), dissolves at 5.0 (edge).
  # Dyad (1,3) forms at 10.0 (inside burn-in), dissolves at 15.0 (inside burn-in).
  # Dyad (1,2) forms at 20.0 (outside), dissolves at 25.0 (outside).
  # Dyad (1,2) forms at 1.0 (inside burn-in), dissolves at 5.0 (edge).
  # Dyad (1,3) forms at 10.0 (inside burn-in), dissolves at 15.0 (inside burn-in).
  # Dyad (1,2) forms at 20.0 (outside), dissolves at 25.0 (outside).
  # Dyad (1,3) forms at 31.0 (outside). 
  events <- matrix(c(
    1.0,  1, 2, 1,
    5.0,  1, 2, 0,
    10.0, 1, 3, 1,
    15.0, 1, 3, 0,
    20.0, 1, 2, 1,
    25.0, 1, 2, 0,
    30.0, 2, 3, 1, # Noise interaction to force a snapshot at t=30
    31.0, 1, 3, 1,
    100.0, 2, 3, 0,
    100.0, 1, 3, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")
  
  # We test 'inertia' because its tracking accumulates historically and persistently validates
  f_0_1 <- ~ inertia
  
  fit1 <- dem(
    events = events,
    formula_0_1 = f_0_1,
    formula_1_0 = ~ Intercept,
    n_nodes = n_nodes,
    directed = TRUE,
    control = control.redeem(build_time = 18.0, return_data = TRUE)
  )
  
  dt <- as.data.table(fit1$model_0_1$data)
  dt_1_0 <- as.data.table(fit1$model_1_0$data)
  
  # Expect that no rows exist with time < 0.0 (relative to build_time=18.0)
  expect_true(all(dt$time >= 0.0))
  
  # Dyad (1,2) is focal in model_1_0 for the interval [20.0, 25.0] where it dissolves.
  dyad_1_2_t25 <- dt_1_0[from == 1 & to == 2 & time <= 25.0 & time_new > 24.99]
  expect_true(nrow(dyad_1_2_t25) > 0)
  # It had one previous formation inside burn-in (t=1.0) and one formation at t=20.0.
  # Snapshooting the dissolution at t=25.0, it has collaborated twice.
  expect_equal(as.numeric(dyad_1_2_t25$inertia), 1)
  
  # Dyad (1,3) is focal in model_1_0 for the interval [30.0, 31.0] where it dissolves? No, it forms at 31.
  # Let's check when (1,3) is focal. 
  # In model_1_0, we have dyad (2,3) focal at 31.0 (noise interaction).
  dyad_2_3_t31 <- dt_1_0[from == 2 & to == 3 & time <= 31.0 & time_new > 30.99]
  expect_true(nrow(dyad_2_3_t31) > 0)
  # inertia for (1,3) in the same snapshot at t=31.0? 
})

