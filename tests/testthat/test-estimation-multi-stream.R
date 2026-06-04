library(testthat)
library(redeem)
library(data.table)

test_that("Multi-stream preprocessing works correctly", {
  # 1. Create two small synthetic event streams
  events_A <- matrix(c(
    1, 1, 5, 1,
    2, 2, 6, 1,
    5, 1, 5, 0,
    6, 2, 6, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events_A) <- c("time", "from", "to", "type")

  events_B <- matrix(c(
    0.5, 1, 5, 1,
    1.5, 1, 5, 0,
    1.5, 2, 6, 1,
    2.5, 2, 6, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events_B) <- c("time", "from", "to", "type")

  n_nodes <- 10

  # 2. Test formula_preprocess
  pre <- redeem:::formula_preprocess(
    formula_0_1 = ~ inertia(event_stream = events_B),
    formula_1_0 = ~Intercept,
    events = events_A,
    n_nodes = n_nodes
  )

  # 3. Test preprocess_multi_stream
  dt_mat <- redeem:::preprocess_multi_stream(
    preprocessed = pre,
    n_nodes = n_nodes,
    verbose = FALSE,
    directed = FALSE,
    simultaneous_interactions = TRUE,
    build_time = 0,
    model_type = "dem"
  )

  dt <- as.data.table(dt_mat)
  # Internal name for inertia (unsuffixed after mapping back)
  expect_true("inertia_identity" %in% names(dt))

  # 4. Test full estimation path
  fit <- dem(
    events = events_A,
    formula_0_1 = ~ inertia(event_stream = events_B),
    formula_1_0 = ~Intercept,
    n_nodes = n_nodes,
    control = control.redeem(estimate = "NR")
  )

  expect_s3_class(fit, "dem")
  # Check if any coefficient starts with inertia
  coef_names <- names(fit$model_0_1$coefficients)
  expect_true(any(grepl("^inertia", coef_names)))
})

test_that("Multi-stream security checks work", {
  events_A <- matrix(c(1, 1, 5, 1, 5, 1, 5, 0), ncol = 4, byrow = TRUE)
  colnames(events_A) <- c("time", "from", "to", "type")
  n_nodes <- 10

  # 1. Node index too high
  events_B_bad_node <- matrix(c(0.5, 1, 11, 1, 1.5, 1, 11, 0), ncol = 4, byrow = TRUE)
  colnames(events_B_bad_node) <- c("time", "from", "to", "type")

  expect_error(
    dem(events = events_A, formula_0_1 = ~ inertia(event_stream = events_B_bad_node), formula_1_0 = ~1, n_nodes = n_nodes),
    "node indices outside the range"
  )

  # 2. Too few columns
  events_B_bad_cols <- matrix(c(0.5, 1, 5), ncol = 3, byrow = TRUE)
  expect_error(
    dem(events = events_A, formula_0_1 = ~ inertia(event_stream = events_B_bad_cols), formula_1_0 = ~1, n_nodes = n_nodes),
    "at least 4 columns"
  )
})

test_that("Multi-stream simulation is prevented", {
  events_A <- matrix(c(1, 1, 5, 1, 5, 1, 5, 0), ncol = 4, byrow = TRUE)
  colnames(events_A) <- c("time", "from", "to", "type")
  events_B <- matrix(c(0.5, 1, 5, 1, 1.5, 1, 5, 0), ncol = 4, byrow = TRUE)
  colnames(events_B) <- c("time", "from", "to", "type")
  n_nodes <- 10

  # Attempt to simulate with multi-stream formula
  expect_error(
    dem.simulate(
      n_nodes = n_nodes,
      formula_0_1 = ~ inertia(event_stream = events_B),
      formula_1_0 = ~1,
      coef_0_1 = 0,
      coef_1_0 = 0,
      n_events = 10
    ),
    "Multi-stream simulation is not yet supported."
  )
})
