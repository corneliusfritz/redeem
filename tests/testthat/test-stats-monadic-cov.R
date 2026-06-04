library(testthat)
library(redeem)

test_that("monadic_cov correctly transforms vector to dyadic matrix", {
  n <- 5
  v <- 1:5
  fun <- function(a, b) a + b

  # Undirected case
  # For undirected, we expect a symmetric matrix where X[i,j] = fun(v[i], v[j])
  # Our implementation applies to upper triangle and mirrors.
  # redeem:::InitRedeemTerm calls are internal, so we test through formula_preprocess or by manually calling the Init logic if we can.
  # Actually, the best way is to test the InitRedeemT- [x] Verify `monadic_cov` with new tests
  # Actually, the best way is to test the redeem:::InitRedeemTerm.monadic_cov directly if it's exported or by calling redeem:::InitRedeemTerm.
  arglist <- list(data = v, fun = fun, base_name = "monadic_cov")
  attr(arglist, "process") <- "0-1"

  # directed = FALSE
  res_undirected <- redeem:::InitRedeemTerm("monadic_cov", arglist, model_type = "rem", process = "0-1", n_nodes = n, directed = FALSE)
  expect_equal(res_undirected$base_name, "monadic_cov")
  mat_undirected <- res_undirected$data
  expect_true(is.matrix(mat_undirected))
  expect_equal(dim(mat_undirected), c(n, n))
  expect_equal(mat_undirected[1, 2], 1 + 2)
  expect_equal(mat_undirected[2, 1], 1 + 2) # Symmetric
  expect_true(isSymmetric(mat_undirected))

  # directed = TRUE
  res_directed <- redeem:::InitRedeemTerm("monadic_cov", arglist, model_type = "rem", process = "0-1", n_nodes = n, directed = TRUE)
  mat_directed <- res_directed$data
  expect_equal(mat_directed[1, 2], 3)
  expect_equal(mat_directed[2, 1], 3)

  # Test with asymmetric function
  fun_asym <- function(a, b) a - b
  arglist_asym <- list(data = v, fun = fun_asym, base_name = "monadic_cov")
  attr(arglist_asym, "process") <- "0-1"

  # directed = TRUE (asymmetric)
  res_asym_dir <- redeem:::InitRedeemTerm("monadic_cov", arglist_asym, model_type = "rem", process = "0-1", n_nodes = n, directed = TRUE)
  expect_equal(res_asym_dir$data[1, 2], 1 - 2)
  expect_equal(res_asym_dir$data[2, 1], 2 - 1)

  # directed = FALSE (symmetric via mirroring upper triangle)
  res_asym_undir <- redeem:::InitRedeemTerm("monadic_cov", arglist_asym, model_type = "rem", process = "0-1", n_nodes = n, directed = FALSE)
  # i=1, j=2: val = 1-2 = -1. mat[1,2]=mat[2,1]=-1.
  expect_equal(res_asym_undir$data[1, 2], -1)
  expect_equal(res_asym_undir$data[2, 1], -1)
  expect_true(isSymmetric(res_asym_undir$data))
})

test_that("monadic_cov handles time-varying data", {
  n <- 3
  v_list <- list("0" = c(1, 1, 1), "10" = c(2, 2, 2))
  fun <- function(a, b) a * b

  arglist <- list(data = v_list, fun = fun, base_name = "monadic_cov")
  attr(arglist, "process") <- "0-1"

  res <- redeem:::InitRedeemTerm("monadic_cov", arglist, model_type = "rem", process = "0-1", n_nodes = n, directed = FALSE)
  expect_true(is.list(res$data))
  expect_equal(names(res$data), c("0", "10"))
  expect_equal(res$data[["0"]][1, 2], 1)
  expect_equal(res$data[["10"]][1, 2], 4)
})
