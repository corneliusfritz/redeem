library(testthat)
library(redeem)

test_that("Windowed statistics change as expected in simulation and preprocessing (Undirected)", {
  n_nodes <- 5
  window_size <- 1
  total_time <- 20

  # 1. Simulate data with windowed common partner statistic
  formula <- ~ common_partner(history="general", window=window_size)
  coef <- 1.0
  names(coef) <- "common_partner(history=\"general\", window=5)"

  set.seed(42)
  sim_events <- rem.simulate(
    formula = formula,
    coef = coef,
    n_nodes = n_nodes,
    time = total_time,
    seed = 42,
    directed = FALSE
  )

  # 2. Fit the model and return preprocessed data
  fit <- rem(
    events = sim_events,
    n_nodes = n_nodes,
    formula = formula,
    directed = FALSE,
    control = control.redeem(return_data = TRUE)
  )
  prep_data <- fit$model$data

  stat_col <- grep("common_partner", colnames(prep_data), value = TRUE)
  expect_true(length(stat_col) > 0)

  # 3. Verify windowing logic on random snapshots
  set.seed(123)
  for (row_idx in sample(seq_len(nrow(prep_data)), 20)) {
    row <- prep_data[row_idx, ]
    t <- row$time
    i <- row$from
    j <- row$to
    observed_val <- row[[stat_col]]

    # Calculate expected value manually from sim_events at time t
    # prep_data snapshots represent the state AFTER events at time 't'
    history <- sim_events[sim_events[, 1] <= (t + 1e-10), , drop = FALSE]
    recent_history <- history[history[, 1] > (t - window_size + 1e-10), , drop = FALSE]

    # Partners of i (undirected)
    partners_i <- unique(c(recent_history[recent_history[, 2] == i, 3],
                           recent_history[recent_history[, 3] == i, 2]))
    # Partners of j (undirected)
    partners_j <- unique(c(recent_history[recent_history[, 2] == j, 3],
                           recent_history[recent_history[, 3] == j, 2]))

    expected_val <- length(intersect(partners_i, partners_j))

    expect_equal(observed_val, expected_val,
                 info = sprintf("Mismatch (Undirected) at t=%f, pair (%d, %d)", t, i, j))
  }
})


test_that("Windowed statistics change as expected in simulation and preprocessing (Directed OSP)", {
  n_nodes <- 5
  window_size <- 5
  total_time <- 20

  # 1. Simulate data with windowed common partner statistic
  formula <- ~ common_partner(history="general", type="OSP", window=window_size)
  coef <- 1.0
  names(coef) <- "common_partner(history=\"general\", type=\"OSP\", window=5)"

  set.seed(42)
  sim_events <- rem.simulate(
    formula = formula,
    coef = coef,
    n_nodes = n_nodes,
    time = total_time,
    seed = 42,
    directed = TRUE
  )

  # 2. Fit the model and return preprocessed data
  fit <- rem(
    events = sim_events,
    n_nodes = n_nodes,
    formula = formula,
    directed = TRUE,
    control = control.redeem(return_data = TRUE)
  )

  prep_data <- fit$model$data
  stat_col <- grep("common_partner", colnames(prep_data), value = TRUE)
  expect_true(length(stat_col) > 0)

  # 3. Verify windowing logic
  set.seed(123)
  for (row_idx in sample(seq_len(nrow(prep_data)), 20)) {
    row <- prep_data[row_idx, ]
    t <- row$time
    i <- row$from
    j <- row$to
    observed_val <- row[[stat_col]]

    history <- sim_events[sim_events[, 1] <= (t + 1e-10), , drop = FALSE]
    recent_history <- history[history[, 1] > (t - window_size + 1e-10), , drop = FALSE]

    # OSP: i -> k and j -> k
    # Partners of i (outgoing)
    partners_i <- unique(recent_history[recent_history[, 2] == i, 3])
    # Partners of j (outgoing)
    partners_j <- unique(recent_history[recent_history[, 2] == j, 3])

    expected_val <- length(intersect(partners_i, partners_j))

    expect_equal(observed_val, expected_val,
                 info = sprintf("Mismatch (Directed OSP) at t=%f, pair (%d, %d)", t, i, j))
  }
})

test_that("Windowed statistics match normal statistics when window is large", {
  n_nodes <- 5
  total_time <- 10

  set.seed(123)
  events <- data.frame(
    time = sort(runif(20, 0, 10)),
    from = sample(1:n_nodes, 20, replace = TRUE),
    to = sample(1:n_nodes, 20, replace = TRUE),
    type = 1
  )
  # Remove self-loops
  events <- events[events$from != events$to, ]
  events <- as.matrix(events)

  # Formula with window = Inf
  formula_inf <- ~ common_partner(history="general", type="OSP", window=Inf)
  # Formula with window = total_time
  formula_win <- ~ common_partner(history="general", type="OSP", window=total_time)

  fit_inf <- rem(
    events = events,
    n_nodes = n_nodes,
    formula = formula_inf,
    directed = TRUE,
    control = control.redeem(return_data = TRUE)
  )

  fit_win <- rem(
    events = events,
    n_nodes = n_nodes,
    formula = formula_win,
    directed = TRUE,
    control = control.redeem(return_data = TRUE)
  )

  stat_col_inf <- grep("common_partner", colnames(fit_inf$model$data), value = TRUE)
  stat_col_win <- grep("common_partner", colnames(fit_win$model$data), value = TRUE)

  expect_equal(fit_inf$model$data[[stat_col_inf]], fit_win$model$data[[stat_col_win]])
})
