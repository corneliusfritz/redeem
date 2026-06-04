library(testthat)
library(redeem)

test_that("check_matrix correctly identifies invalid data", {
  # 1. Valid sequence
  m1 <- matrix(c(
    1.0, 1, 2, 1,
    2.0, 1, 2, 0
  ), ncol = 4, byrow = TRUE)
  expect_true(check_matrix(m1))

  # 2. Overlap (Still an error)
  m2 <- matrix(c(
    1.0, 1, 2, 1,
    1.5, 1, 2, 1
  ), ncol = 4, byrow = TRUE)
  w2 <- capture_warnings(res2 <- check_matrix(m2))
  expect_false(res2)
  expect_match(paste(w2, collapse = " "), "Overlap found")
  
  # Check return_matrix=TRUE with errors
  w2_mat <- capture_warnings({ res2_mat <- check_matrix(m2, return_matrix = TRUE) })
  expect_true(is.matrix(res2_mat))
  expect_true(!is.null(attr(res2_mat, "has_errors")))

  # 3. End without start (Silent Auto-repaired)
  m3 <- matrix(c(
    1.0, 1, 2, 0
  ), ncol = 4, byrow = TRUE)
  # Should NOT throw a warning anymore as requested by user
  res3 <- expect_silent(check_matrix(m3, return_matrix = TRUE))
  expect_true(is.matrix(res3))
  expect_equal(nrow(res3), 2)
  expect_equal(res3[1, 4], 1) # Added start (at t=1.0)
  expect_equal(res3[2, 4], 0) # Original end

  # 4. Unclosed interaction (Silent keep)
  m4 <- matrix(c(
    1.0, 1, 2, 1
  ), ncol = 4, byrow = TRUE)
  # Should NOT throw a warning anymore as requested by user
  res4 <- expect_silent(check_matrix(m4, return_matrix = TRUE))
  expect_true(is.matrix(res4))
  expect_equal(nrow(res4), 1)
  
  # 5. Missing values (Error)
  m5 <- matrix(c(
    1.0, 1, 2, 1,
    NA, 1, 2, 0
  ), ncol = 4, byrow = TRUE)
  w5 <- capture_warnings(res5 <- check_matrix(m5))
  expect_false(res5)
  expect_match(paste(w5, collapse = " "), "Found missing values")
})

test_that("check_matrix handles character IDs and formatting", {
  m_char <- matrix(c(
    "1.0", "ActorA", "ActorB", "1",
    "2.0", "ActorA", "ActorB", "1"
  ), ncol = 4, byrow = TRUE)
  # Only overlap should warn/error
  w_char <- capture_warnings(res_char <- check_matrix(m_char))
  expect_false(res_char)
  expect_match(paste(w_char, collapse = " "), "Overlap found for dyad \\(ActorA, ActorB\\)")
})

test_that("dem() handles unclosed data gracefully (censored)", {
  m_unclosed <- matrix(c(
    1.0, 1, 2, 1
  ), ncol = 4, byrow = TRUE)
  colnames(m_unclosed) <- c("time", "from", "to", "type")
  
  # Should now work silently
  fit <- expect_silent(dem(events = m_unclosed, n_nodes = 5, formula_0_1 = ~ 1, control = control.redeem(it_max = 1)))
  expect_s3_class(fit, "dem")
})

test_that("dem() still aborts on severe errors (Overlap)", {
  m_bad <- matrix(c(
    1.0, 1, 2, 1,
    1.5, 1, 2, 1
  ), ncol = 4, byrow = TRUE)
  colnames(m_bad) <- c("time", "from", "to", "type")
  
  expect_error(
    suppressWarnings(dem(events = m_bad, n_nodes = 5, formula_0_1 = ~ 1, control = control.redeem(check_matrix = TRUE))),
    "The provided event data contains fatal errors"
  )
})

test_that("check_matrix handles multi-error reporting", {
  m_multi <- matrix(c(
    1.0, 1, 2, 1,
    1.5, 1, 2, 1, # Error 1 (Overlap)
    2.0, 1, 2, 0,
    1.0, 3, 4, 0, # Silently Repaired
    3.0, 5, 6, 1  # Silently Kept
  ), ncol = 4, byrow = TRUE)
  
  w <- capture_warnings(res <- check_matrix(m_multi))
  expect_false(res)
  # Only Overlap should be in the warnings
  expect_match(paste(w, collapse = " "), "Overlap found")
  expect_no_match(paste(w, collapse = " "), "Excluding")
})
