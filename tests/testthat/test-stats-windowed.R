library(testthat)
library(redeem)

test_that("windowed inertia behaves correctly", {
  n <- 3
  # Dyad (1,2) interaction at time 1
  events <- matrix(c(
    1, 1, 2, 1,
    2, 1, 2, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")

  # Window = 5. The interaction at time 1 should "expire" at time 6.
  # We'll check the statistics at time 2 and time 7.
  
  # To check statistics manually, we can use preprocess
  # We need to add an event at time 7 to see the state then.
  events_with_check <- rbind(events, c(7, 1, 3, 1))
  
  # R side formula preprocessing
  pre <- redeem:::formula_preprocess(
    formula_0_1 = ~ inertia(window = 5),
    events = events_with_check,
    n_nodes = n,
    directed = TRUE
  )
  
  # C++ side preprocessing
  processed <- redeem:::preprocess(
    edgelist = as.matrix(pre$events),
    terms = pre$term_names,
    data_list = pre$data_list,
    transformations = pre$transformation_list,
    n_nodes = n,
    verbose = FALSE,
    directed = TRUE,
    simultaneous_interactions = FALSE,
    window_map = pre$window_map
  )
  
  # processed columns: time, from, to, type, status, event, from_avail, to_avail, Intercept, inertia_wt11
  # The column index for inertia_wt11 should be 10 (1-based) if Intercept is 9.
  # Wait, let's find the column.
  topo_names <- c("time_new", "time", "pair_id", "status", "event", "from", "to", "from_avail", "to_avail")
  processed <- as.data.frame(processed)
  colnames(processed) <- c(topo_names, as.character(pre$term_names))
  col_idx <- grep("inertia_wt", colnames(processed))
  
  pair12 <- 1 
  
  # Event 1: time 1, (1,2,1). inertia before = 0
  # Interval [0, 1), time_new == 1
  val1 <- processed[processed$pair_id == pair12 & abs(processed$time_new - 1) < 1e-7, col_idx]
  expect_equal(as.numeric(val1), 0)
  
  # Event 2: time 2, (1,2,0). inertia before = 1
  # Interval [1, 2), time_new == 2
  val2 <- processed[processed$pair_id == pair12 & abs(processed$time_new - 2) < 1e-7, col_idx]
  expect_equal(as.numeric(val2), 1)
  
  # Event 3: time 7, (1,3,1). inertia (1,2) should be 0 because 1+5=6 < 7.
  # Interval [6, 7), time_new == 7
  val3 <- processed[processed$pair_id == pair12 & abs(processed$time_new - 7) < 1e-7, col_idx]
  expect_equal(as.numeric(val3), 0)
  
  # To check inertia for (1,2) at time 5, we need an event at time 5.
  events_check_5 <- rbind(events, c(5, 1, 2, 1))
  pre3 <- redeem:::formula_preprocess(
    formula_0_1 = ~ inertia(window = 5),
    events = events_check_5,
    n_nodes = n,
    directed = TRUE
  )
  processed3 <- redeem:::preprocess(
    edgelist = as.matrix(pre3$events),
    terms = pre3$term_names,
    data_list = pre3$data_list,
    transformations = pre3$transformation_list,
    n_nodes = n,
    verbose = FALSE,
    directed = TRUE,
    simultaneous_interactions = FALSE,
    window_map = pre3$window_map
  )
  processed3 <- as.data.frame(processed3)
  colnames(processed3) <- c(topo_names, as.character(pre3$term_names))
  col_idx3 <- grep("inertia_wt", colnames(processed3))
  # Interval [?, 5), time_new == 5
  val5 <- processed3[processed3$pair_id == pair12 & abs(processed3$time_new - 5) < 1e-7, col_idx3]
  expect_equal(as.numeric(val5), 1)
})

test_that("windowed degree behaves correctly", {
  n <- 3
  # (1,2) at time 1, (1,3) at time 2
  events <- matrix(c(
    1, 1, 2, 1,
    1.5, 1, 2, 0,
    2, 1, 3, 1,
    2.5, 1, 3, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")
  
  topo_names <- c("time_new", "time", "pair_id", "status", "event", "from", "to", "from_avail", "to_avail")
  
  check_times <- c(1, 2, 4, 5)
  results <- sapply(check_times, function(tt) {
    ev <- rbind(events, c(tt, 1, 2, 1))
    pre <- redeem:::formula_preprocess(~ degree(type="out_sender", window=2), events=ev, n_nodes=n, directed=TRUE)
    proc <- redeem:::preprocess(as.matrix(pre$events), pre$term_names, pre$data_list, pre$transformation_list, n, FALSE, TRUE, FALSE, pre$window_map)
    proc <- as.data.frame(proc)
    colnames(proc) <- c(topo_names, as.character(pre$term_names))
    col <- grep("degree_out_sender_wt", colnames(proc))
    # Filter by time_new == tt
    as.numeric(proc[which.min(abs(proc$time_new - tt)), col])
  })
  
  expect_equal(as.numeric(results), c(0, 1, 1, 0))
})

test_that("manual calculation of windowed general_common_partner is correct", {
  n <- 4
  # Scenario for OSP (Outward Shared Partner): 1->k and 2->k
  # t=1: 1->3
  # t=2: 2->3
  # t=3: 1->4
  # t=4: 2->4
  # t=4.5: dummy event (1->2) to check state. common_partner of (1,2) should be 2.
  # t=6.5: dummy event (1->2). 1->3 expired at 6. common_partner of (1,2) should be 1.
  # t=7.5: dummy event (1->2). 2->3 expired at 7. common_partner of (1,2) should be 1 (because 1->4 and 2->4 are still active).
  # t=8.5: dummy event (1->2). 1->4 expired at 8. common_partner of (1,2) should be 0.
  
  events <- matrix(c(
    1, 1, 3, 1,
    2, 2, 3, 1,
    3, 1, 4, 1,
    4, 2, 4, 1,
    4.5, 1, 2, 1,
    6.5, 1, 2, 1,
    7.5, 1, 2, 1,
    8.5, 1, 2, 1
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")
  
  pre <- redeem:::formula_preprocess(
    formula_0_1 = ~ general_common_partner(window = 5),
    events = events,
    n_nodes = n,
    directed = TRUE
  )
  
  processed <- redeem:::preprocess(
    edgelist = as.matrix(pre$events),
    terms = pre$term_names,
    data_list = pre$data_list,
    transformations = pre$transformation_list,
    n_nodes = n,
    verbose = FALSE,
    directed = TRUE,
    simultaneous_interactions = FALSE,
    window_map = pre$window_map
  )
  
  topo_names <- c("time_new", "time", "pair_id", "status", "event", "from", "to", "from_avail", "to_avail")
  processed <- as.data.frame(processed)
  colnames(processed) <- c(topo_names, as.character(pre$term_names))
  
  col_idx <- grep("general_common_partner", colnames(processed))
  
  # Expected values at the ends of intervals for from=1, to=2:
  # The row for interval ending at 4.5
  val_4_5 <- processed[processed$from == 1 & processed$to == 2 & abs(processed$time_new - 4.5) < 1e-7, col_idx]
  expect_equal(as.numeric(val_4_5), 2)
  
  # The row for interval ending at 6.5
  val_6_5 <- processed[processed$from == 1 & processed$to == 2 & abs(processed$time_new - 6.5) < 1e-7, col_idx]
  expect_equal(as.numeric(val_6_5), 1)
  
  # The row for interval ending at 7.5
  val_7_5 <- processed[processed$from == 1 & processed$to == 2 & abs(processed$time_new - 7.5) < 1e-7, col_idx]
  expect_equal(as.numeric(val_7_5), 1)
  
  # The row for interval ending at 8.5
  val_8_5 <- processed[processed$from == 1 & processed$to == 2 & abs(processed$time_new - 8.5) < 1e-7, col_idx]
  expect_equal(as.numeric(val_8_5), 0)
})

test_that("windowed general_count_absdiff behaves correctly under event expiration", {
  n <- 4
  # t=1: 1->3
  # t=2: 2->3
  # t=3: 2->3
  # t=4.5: 1->2 (dummy event to check)
  # t=6.5: 1->2 (dummy event to check)
  # t=7.5: 1->2 (dummy event to check)
  # t=8.5: 1->2 (dummy event to check)
  events <- matrix(c(
    1.0, 1, 3, 1,
    2.0, 2, 3, 1,
    3.0, 2, 3, 1,
    4.5, 1, 2, 1,
    6.5, 1, 2, 1,
    7.5, 1, 2, 1,
    8.5, 1, 2, 1
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")

  pre <- redeem:::formula_preprocess(
    formula_0_1 = ~ general_count_absdiff(window = 5),
    events = events,
    n_nodes = n,
    directed = FALSE
  )

  processed <- redeem:::preprocess(
    edgelist = as.matrix(pre$events),
    terms = pre$term_names,
    data_list = pre$data_list,
    transformations = pre$transformation_list,
    n_nodes = n,
    verbose = FALSE,
    directed = FALSE,
    simultaneous_interactions = FALSE,
    window_map = pre$window_map
  )

  topo_names <- c("time_new", "time", "pair_id", "status", "event", "from", "to", "from_avail", "to_avail")
  processed <- as.data.frame(processed)
  colnames(processed) <- c(topo_names, as.character(pre$term_names))

  col_idx <- grep("general_count_absdiff_wt", colnames(processed))

  # t=4.5 check: out_count(1)=1, out_count(2)=2. absdiff = 1.
  val_4_5 <- processed[processed$from == 1 & processed$to == 2 & abs(processed$time_new - 4.5) < 1e-7, col_idx]
  expect_equal(as.numeric(val_4_5), 1)

  # t=6.5 check: 1->3 expired at 6.0. out_count(1)=0, out_count(2)=2. absdiff = 2.
  val_6_5 <- processed[processed$from == 1 & processed$to == 2 & abs(processed$time_new - 6.5) < 1e-7, col_idx]
  expect_equal(as.numeric(val_6_5), 2)

  # t=7.5 check: 2->3 at t=2.0 expired at 7.0. out_count(1)=0, out_count(2)=1. absdiff = 1.
  val_7_5 <- processed[processed$from == 1 & processed$to == 2 & abs(processed$time_new - 7.5) < 1e-7, col_idx]
  expect_equal(as.numeric(val_7_5), 1)

  # t=8.5 check: 2->3 at t=3.0 expired at 8.0. out_count(1)=0, out_count(2)=0. absdiff = 0.
  val_8_5 <- processed[processed$from == 1 & processed$to == 2 & abs(processed$time_new - 8.5) < 1e-7, col_idx]
  expect_equal(as.numeric(val_8_5), 0)
})
