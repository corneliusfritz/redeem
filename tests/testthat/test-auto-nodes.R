library(testthat)
library(redeem)

test_that("process_event_actors helper works correctly", {
  # 1. Numeric actors, n_nodes = NULL, directed = TRUE
  events_num <- data.frame(time = 1:3, from = c(2, 5, 3), to = c(4, 1, 2), type = 1)
  res <- redeem:::process_event_actors(events_num, n_nodes = NULL, directed = TRUE)
  expect_equal(res$n_nodes, 5)
  expect_equal(res$events$from, c(2, 5, 3))
  expect_equal(res$events$to, c(4, 1, 2))

  # 2. String actors, n_nodes = NULL, directed = TRUE
  events_str <- data.frame(time = 1:3, from = c("Bob", "Charlie", "Bob"), to = c("Alice", "Bob", "Charlie"), type = 1)
  res <- redeem:::process_event_actors(events_str, n_nodes = NULL, directed = TRUE)
  # Unique sorted: Alice, Bob, Charlie (1, 2, 3)
  expect_equal(res$n_nodes, 3)
  expect_equal(res$events$from, c(2, 3, 2)) # Bob=2, Charlie=3, Bob=2
  expect_equal(res$events$to, c(1, 2, 3))   # Alice=1, Bob=2, Charlie=3

  # 3. Numeric actors, undirected (directed = FALSE)
  events_num_undir <- data.frame(time = 1:2, from = c(5, 1), to = c(2, 4), type = 1)
  res <- redeem:::process_event_actors(events_num_undir, n_nodes = NULL, directed = FALSE)
  expect_equal(res$n_nodes, 5)
  # 5 > 2 -> swap to 2, 5; 1 < 4 -> keep 1, 4
  expect_equal(res$events$from, c(2, 1))
  expect_equal(res$events$to, c(5, 4))

  # 4. String actors, undirected (directed = FALSE)
  # Unique: Alice (1), Bob (2)
  # Row 1: Bob (2) -> Alice (1). Since undirected and 2 > 1, swap to Alice (1) -> Bob (2)
  events_str_undir <- data.frame(time = 1, from = "Bob", to = "Alice", type = 1)
  res <- redeem:::process_event_actors(events_str_undir, n_nodes = NULL, directed = FALSE)
  expect_equal(res$n_nodes, 2)
  expect_equal(res$events$from, 1)
  expect_equal(res$events$to, 2)

  # 5. Stop if more actors than provided n_nodes
  expect_error(
    redeem:::process_event_actors(events_num, n_nodes = 3, directed = TRUE),
    "More unique actors.*identified.*provided in the n_nodes"
  )
  expect_error(
    redeem:::process_event_actors(events_str, n_nodes = 2, directed = TRUE),
    "More unique actors.*identified.*provided in the n_nodes"
  )

  # 6. Stop if max actor ID is greater than n_nodes (for numeric)
  # Unique count is 3 (1, 2, 5), which is <= n_nodes (4), but max actor ID (5) > 4
  events_num_gap <- data.frame(time = 1:2, from = c(1, 2), to = c(5, 1), type = 1)
  expect_error(
    redeem:::process_event_actors(events_num_gap, n_nodes = 4, directed = TRUE),
    "Maximum actor ID.*larger than the provided n_nodes"
  )
})

test_that("dem and rem work with automatic n_nodes and string actors", {
  events_str <- matrix(c(
    1.2, "Alice", "Bob", 1,
    2.5, "Alice", "Bob", 0,
    3.1, "Bob", "Charlie", 1,
    4.4, "Bob", "Charlie", 0
  ), ncol = 4, byrow = TRUE)
  colnames(events_str) <- c("time", "from", "to", "type")

  # 1. dem with character actors and n_nodes = NULL
  fit_dem <- dem(
    events = events_str,
    formula_0_1 = ~1,
    formula_1_0 = ~1,
    n_nodes = NULL
  )
  expect_s3_class(fit_dem, "dem")
  expect_equal(fit_dem$n_nodes, 3)

  # 2. rem with character actors and n_nodes = NULL
  fit_rem <- rem(
    events = events_str,
    formula = ~1,
    n_nodes = NULL
  )
  expect_s3_class(fit_rem, "rem")
  expect_equal(fit_rem$n_nodes, 3)

  # 3. Numeric events with n_nodes = NULL
  events_num <- matrix(c(
    1.2, 1, 3, 1,
    2.5, 1, 3, 0,
    3.1, 2, 4, 1,
    4.4, 2, 4, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events_num) <- c("time", "from", "to", "type")

  fit_dem_num <- dem(
    events = events_num,
    formula_0_1 = ~1,
    formula_1_0 = ~1,
    n_nodes = NULL
  )
  expect_equal(fit_dem_num$n_nodes, 4)

  fit_rem_num <- rem(
    events = events_num,
    formula = ~1,
    n_nodes = NULL
  )
  expect_equal(fit_rem_num$n_nodes, 4)
})
