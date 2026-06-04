library(testthat)
library(redeem)

test_that("dem.simulate and dem: Only baseline", {
  set.seed(123)
  n_nodes <- 10
  cp_seq <- c(100, 200, 300)
  # 4 slices (0-100, 100-200, 200-300, 300+)
  true_baseline_0_1 <- c(-0.5, 0, 0.5)
  true_baseline_1_0 <- c(-1, -0.5, 0)

  # Simulate (Intercept is automatically added, auto-shifting will handle full baseline)
  events <- dem.simulate(
    formula_0_1 = ~ baseline(changepoints = cp_seq),
    formula_1_0 = ~ baseline(changepoints = cp_seq),
    baseline_0_1 = true_baseline_0_1,
    baseline_1_0 = true_baseline_1_0,
    time = 400,
    n_nodes = n_nodes,
    directed = FALSE
  )
  # Estimate
  fit <- dem(
    events = events,
    formula_0_1 = ~ baseline(changepoints = cp_seq),
    formula_1_0 = ~ baseline(changepoints = cp_seq),
    n_nodes = n_nodes,
    directed = FALSE
  )

  # Check recovery (absolute baseline levels when no degrees/intercept)
  expect_equal(as.vector(fit$model_0_1$est_time), as.vector(c(0, true_baseline_0_1)), tolerance = 0.5)
  expect_equal(as.vector(fit$model_1_0$est_time), as.vector(c(0, true_baseline_1_0)), tolerance = 0.5)
})

test_that("dem.simulate and dem: Only degrees", {
  set.seed(123)
  n_nodes <- 10
  true_degree <- rnorm(n_nodes, 0, 0.5)
  true_degree <- true_degree - mean(true_degree)

  # Simulate
  events <- dem.simulate(
    formula_0_1 = ~degrees,
    formula_1_0 = ~degrees,
    coef_degree_0_1 = true_degree,
    coef_degree_1_0 = true_degree,
    time = 400,
    n_nodes = n_nodes,
    directed = FALSE
  )
  # Estimate
  fit <- dem(
    events = events,
    formula_0_1 = ~degrees,
    formula_1_0 = ~degrees,
    n_nodes = n_nodes,
    directed = FALSE
  )

  expect_gt(cor(fit$model_0_1$est_degree, true_degree), 0.8)
  expect_gt(cor(fit$model_1_0$est_degree, true_degree), 0.8)
})

test_that("dem.simulate and dem: Baseline and degrees", {
  set.seed(123)
  n_nodes <- 10
  cp_seq <- c(100, 200, 300)
  true_baseline_0_1 <- c(0.5, 1.0, 1.5)
  true_baseline_1_0 <- c(0.5, 1.0, 1.5)
  true_degree <- rnorm(n_nodes, 0, 0.5)

  # Simulate (degrees present, baseline must be offsets: length 3)
  events <- dem.simulate(
    formula_0_1 = ~ baseline(changepoints = cp_seq) + degrees,
    formula_1_0 = ~ baseline(changepoints = cp_seq) + degrees,
    baseline_0_1 = true_baseline_0_1,
    baseline_1_0 = true_baseline_1_0,
    coef_degree_0_1 = true_degree,
    coef_degree_1_0 = true_degree,
    time = 600,
    n_nodes = n_nodes,
    directed = FALSE
  )

  # Estimate
  fit <- dem(
    events = events,
    formula_0_1 = ~ baseline(changepoints = cp_seq) + degrees,
    formula_1_0 = ~ baseline(changepoints = cp_seq) + degrees,
    n_nodes = n_nodes,
    directed = FALSE
  )

  expect_gt(cor(fit$model_0_1$est_degree, true_degree), 0.7)
  # Check recovery
  expect_equal(as.vector(fit$model_0_1$est_time), as.vector(true_baseline_0_1), tolerance = 0.5)
})

test_that("rem.simulate and rem: Only baseline", {
  set.seed(123)
  n_nodes <- 10
  cp_seq <- c(100, 200)
  true_baseline <- c(-1, 0)

  # Simulate (Intercept automatically added, auto-shift handles full baseline)
  events <- rem.simulate(
    formula = ~ baseline(changepoints = cp_seq),
    baseline = true_baseline,
    time = 400,
    n_nodes = n_nodes,
    directed = FALSE
  )
  # Estimate
  fit <- rem(
    events = events,
    formula = ~ baseline(changepoints = cp_seq),
    n_nodes = n_nodes,
    directed = FALSE
  )

  # Check recovery (baseline differences)
  est_baseline <- if (fit$model$full_baseline) fit$model$est_time else c(0, fit$model$est_time)
  expect_equal(as.vector(diff(est_baseline)), as.vector(diff(c(0, true_baseline))), tolerance = 0.5)
})


test_that("rem.simulate and rem: Baseline and degrees", {
  set.seed(123)
  n_nodes <- 10
  cp_seq <- c(100, 200, 300)
  true_baseline <- c(0.5, 1, 1.5)

  # Simulate (degrees present, baseline must be offsets: length 3)
  events <- rem.simulate(
    formula = ~ degrees + baseline(changepoints = cp_seq),
    baseline = true_baseline,
    time = 400,
    n_nodes = n_nodes,
    directed = FALSE
  )
  # Estimate
  fit <- rem(
    events = events,
    formula = ~ degrees + baseline(changepoints = cp_seq),
    n_nodes = n_nodes,
    directed = FALSE
  )

  # Check recovery
  expect_equal(as.vector(fit$model$est_time), as.vector(true_baseline), tolerance = 0.5)
})
