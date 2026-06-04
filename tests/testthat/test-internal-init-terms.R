test_that("redeem:::InitRedeemTerm dispatcher and validator coverage", {
  # Valid dispatch
  arglist <- list(base_name = "Intercept")
  attr(arglist, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm("Intercept", arglist, model_type = "dem", process = "0-1", n_nodes = 5), "list")
  
  # Invalid term
  expect_error(redeem:::InitRedeemTerm("NonExistentTerm", arglist, "dem", "0-1", 5), "not recognized")
  
  # redeem:::check.RedeemTerm - invalid process
  arglist_err <- list(base_name = "current_interaction")
  attr(arglist_err, "process") <- "0-1"
  expect_error(redeem:::check.RedeemTerm(arglist_err, allowed_processes = "1-0"), "not allowed")
  
  # redeem:::check.RedeemTerm - defaults and types
  arglist_types <- list(base_name = "test", n = "not_numeric", m = 1:5, l = list(a=1))
  attr(arglist_types, "process") <- "0-1"
  
  expect_error(redeem:::check.RedeemTerm(arglist_types, expected = list(n = "numeric")), "must be numeric")
  arglist_types$m <- "not_a_matrix"
  expect_error(redeem:::check.RedeemTerm(arglist_types, expected = list(m = "matrix")), "must be a matrix")
  
  # redeem:::check.RedeemTerm - character spec
  arg_char <- list(base_name="t", val="C")
  attr(arg_char, "process") <- "0-1"
  expect_error(redeem:::check.RedeemTerm(arg_char, expected=list(val=c("A", "B"))), "must be one of")
  
  # matrix_or_list, numeric_or_list
  arg_mat <- list(base_name="t", val="string")
  attr(arg_mat, "process") <- "0-1"
  expect_error(redeem:::check.RedeemTerm(arg_mat, expected=list(val="matrix_or_list")), "matrix, numeric vector, or list")
  expect_error(redeem:::check.RedeemTerm(arg_mat, expected=list(val="numeric_or_list")), "numeric vector or list")
})

test_that("Standard term initializers coverage", {
  # Intercept
  arg_int <- list(base_name="Intercept")
  attr(arg_int, "process") <- "0-1"
  expect_equal(redeem:::InitRedeemTerm.Intercept(arg_int, 5, model_type = "dem", directed = TRUE)$base_name, "Intercept")
  
  # Inertia / reciprocity
  arg_ine <- list(base_name="inertia")
  attr(arg_ine, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm.inertia(arg_ine, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.reciprocity(arg_ine, 5, model_type = "dem", directed = TRUE), "list")
  
  # current_interaction (requires 1-0)
  arg_curr <- list(base_name="current_interaction")
  attr(arg_curr, "process") <- "1-0"
  expect_type(redeem:::InitRedeemTerm.current_interaction(arg_curr, 5, model_type = "dem", directed = TRUE), "list")
  
  # common_partner
  arg_cp <- list(base_name="gcp")
  attr(arg_cp, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm.general_common_partner(arg_cp, 5, model_type = "dem", directed = FALSE), "list")
  expect_type(redeem:::InitRedeemTerm.general_common_partner(arg_cp, 5, model_type = "dem", directed = TRUE), "list")
  
  expect_type(redeem:::InitRedeemTerm.current_common_partner(arg_cp, 5, model_type = "dem", directed = FALSE), "list")
  expect_type(redeem:::InitRedeemTerm.current_common_partner(arg_cp, 5, model_type = "dem", directed = TRUE), "list")
})

test_that("Degree and triangle dispatch coverage", {
  # degree
  arg_deg <- list(history = "general", type = "out_sender")
  attr(arg_deg, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm.degree(arg_deg, n_nodes = 5, model_type = "dem", directed = TRUE), "list")
  
  # count
  expect_type(redeem:::InitRedeemTerm.count(arg_deg, n_nodes = 5, model_type = "dem", directed = TRUE), "list")
  
  # Invalid dispatch
  expect_error(redeem:::InitRedeemTerm.degree(list(history = "invalid"), n_nodes = 5, model_type = "dem", directed = TRUE), "not found")
  
  # triangle
  arg_tri <- list(history="general")
  attr(arg_tri, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm.triangle(arg_tri, n_nodes = 5, model_type = "dem", directed = TRUE), "list")
  expect_error(redeem:::InitRedeemTerm.triangle(list(history="invalid", process="0-1"), n_nodes = 5, model_type = "dem", directed = TRUE), "not found")
  
  # common_partner dispatch
  expect_type(redeem:::InitRedeemTerm.common_partner(arg_tri, n_nodes = 5, model_type = "dem", directed = TRUE), "list")
  
  # triangle specific
  arg_tri_spec <- list(base_name="gt")
  attr(arg_tri_spec, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm.general_triangle(arg_tri_spec, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.current_triangle(arg_tri_spec, 5, model_type = "dem", directed = TRUE), "list")
})

test_that("Directed/Undirected constraints coverage", {
  arg_base <- list()
  attr(arg_base, "process") <- "0-1"
  
  # Directed terms on undirected model
  expect_error(redeem:::InitRedeemTerm.degree_out_sender(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.degree_out_receiver(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.degree_in_sender(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.degree_in_receiver(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  
  expect_error(redeem:::InitRedeemTerm.general_count_out_sender(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.general_count_out_receiver(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.general_count_in_sender(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.general_count_in_receiver(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")

  expect_error(redeem:::InitRedeemTerm.general_degree_out_sender(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.general_degree_out_receiver(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.general_degree_in_sender(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.general_degree_in_receiver(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")

  expect_error(redeem:::InitRedeemTerm.current_degree_out_sender(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.current_degree_out_receiver(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.current_degree_in_sender(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.current_degree_in_receiver(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")

  expect_error(redeem:::InitRedeemTerm.current_count_out_sender(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.current_count_out_receiver(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.current_count_in_sender(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")
  expect_error(redeem:::InitRedeemTerm.current_count_in_receiver(arg_base, 5, model_type = "dem", directed = FALSE), "only available for directed")

  # Undirected terms on directed model
  expect_error(redeem:::InitRedeemTerm.degree_sum(arg_base, 5, model_type = "dem", directed = TRUE), "only available for undirected")
  expect_error(redeem:::InitRedeemTerm.degree_absdiff(arg_base, 5, model_type = "dem", directed = TRUE), "only available for undirected")
  expect_error(redeem:::InitRedeemTerm.general_count_sum(arg_base, 5, model_type = "dem", directed = TRUE), "only available for undirected")
  expect_error(redeem:::InitRedeemTerm.general_count_absdiff(arg_base, 5, model_type = "dem", directed = TRUE), "only available for undirected")
  expect_error(redeem:::InitRedeemTerm.general_degree_sum(arg_base, 5, model_type = "dem", directed = TRUE), "only available for undirected")
  expect_error(redeem:::InitRedeemTerm.general_degree_absdiff(arg_base, 5, model_type = "dem", directed = TRUE), "only available for undirected")
  expect_error(redeem:::InitRedeemTerm.current_degree_sum(arg_base, 5, model_type = "dem", directed = TRUE), "only available for undirected")
  expect_error(redeem:::InitRedeemTerm.current_degree_absdiff(arg_base, 5, model_type = "dem", directed = TRUE), "only available for undirected")
  expect_error(redeem:::InitRedeemTerm.current_count_sum(arg_base, 5, model_type = "dem", directed = TRUE), "only available for undirected")
  expect_error(redeem:::InitRedeemTerm.current_count_absdiff(arg_base, 5, model_type = "dem", directed = TRUE), "only available for undirected")
})

test_that("Directed degree/count terms success paths", {
  arg_base <- list()
  attr(arg_base, "process") <- "0-1"
  
  # Degree terms
  expect_type(redeem:::InitRedeemTerm.degree_out_sender(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.degree_out_receiver(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.degree_in_sender(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.degree_in_receiver(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  
  # General count terms
  expect_type(redeem:::InitRedeemTerm.general_count_out_sender(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.general_count_out_receiver(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.general_count_in_sender(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.general_count_in_receiver(arg_base, 5, model_type = "dem", directed = TRUE), "list")

  # General degree terms
  expect_type(redeem:::InitRedeemTerm.general_degree_out_sender(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.general_degree_out_receiver(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.general_degree_in_sender(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.general_degree_in_receiver(arg_base, 5, model_type = "dem", directed = TRUE), "list")

  # Current degree terms
  expect_type(redeem:::InitRedeemTerm.current_degree_out_sender(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.current_degree_out_receiver(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.current_degree_in_sender(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.current_degree_in_receiver(arg_base, 5, model_type = "dem", directed = TRUE), "list")

  # Current count terms
  expect_type(redeem:::InitRedeemTerm.current_count_out_sender(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.current_count_out_receiver(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.current_count_in_sender(arg_base, 5, model_type = "dem", directed = TRUE), "list")
  expect_type(redeem:::InitRedeemTerm.current_count_in_receiver(arg_base, 5, model_type = "dem", directed = TRUE), "list")
})

test_that("Undirected degree/count terms success paths", {
  arg_base <- list()
  attr(arg_base, "process") <- "0-1"

  expect_type(redeem:::InitRedeemTerm.degree_sum(arg_base, 5, model_type = "dem", directed = FALSE), "list")
  expect_type(redeem:::InitRedeemTerm.degree_absdiff(arg_base, 5, model_type = "dem", directed = FALSE), "list")
  expect_type(redeem:::InitRedeemTerm.general_count_sum(arg_base, 5, model_type = "dem", directed = FALSE), "list")
  expect_type(redeem:::InitRedeemTerm.general_count_absdiff(arg_base, 5, model_type = "dem", directed = FALSE), "list")
  expect_type(redeem:::InitRedeemTerm.general_degree_sum(arg_base, 5, model_type = "dem", directed = FALSE), "list")
  expect_type(redeem:::InitRedeemTerm.general_degree_absdiff(arg_base, 5, model_type = "dem", directed = FALSE), "list")
  expect_type(redeem:::InitRedeemTerm.current_degree_sum(arg_base, 5, model_type = "dem", directed = FALSE), "list")
  expect_type(redeem:::InitRedeemTerm.current_degree_absdiff(arg_base, 5, model_type = "dem", directed = FALSE), "list")
  expect_type(redeem:::InitRedeemTerm.current_count_sum(arg_base, 5, model_type = "dem", directed = FALSE), "list")
  expect_type(redeem:::InitRedeemTerm.current_count_absdiff(arg_base, 5, model_type = "dem", directed = FALSE), "list")
})


test_that("dyadic_cov coverage", {
  arg_base <- list()
  attr(arg_base, "process") <- "0-1"

  # Missing data
  expect_error(redeem:::InitRedeemTerm.dyadic_cov(arg_base, 5, model_type = "dem", directed = TRUE), "requires a 'data' argument")
  
  # Scalar data
  arg_sc <- list(data = 1)
  attr(arg_sc, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm.dyadic_cov(arg_sc, 5, model_type = "dem", directed = TRUE), "list")
  
  # Matrix data
  mat <- matrix(0, 5, 5)
  arg_mat <- list(data = mat)
  attr(arg_mat, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm.dyadic_cov(arg_mat, 5, model_type = "dem", directed = TRUE), "list")
  
  arg_bad_mat <- list(data = matrix(0, 4, 4))
  attr(arg_bad_mat, "process") <- "0-1"
  expect_error(redeem:::InitRedeemTerm.dyadic_cov(arg_bad_mat, 5, model_type = "dem", directed = TRUE), "must be a 5 x 5 matrix")
  
  arg_vec <- list(data = 1:2)
  attr(arg_vec, "process") <- "0-1"
  expect_error(redeem:::InitRedeemTerm.dyadic_cov(arg_vec, 5, model_type = "dem", directed = TRUE), "vector of length 2")

  # List data
  l_data <- list("0" = mat, "10" = mat)
  arg_list <- list(data = l_data)
  attr(arg_list, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm.dyadic_cov(arg_list, 5, model_type = "dem", directed = TRUE), "list")
  
  # List data with change_points
  arg_cp <- list(data = list(mat, mat), change_points = c(0, 10))
  attr(arg_cp, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm.dyadic_cov(arg_cp, 5, model_type = "dem", directed = TRUE), "list")
  
  # List data errors
  arg_no_name <- list(data = list(mat, mat))
  attr(arg_no_name, "process") <- "0-1"
  expect_error(redeem:::InitRedeemTerm.dyadic_cov(arg_no_name, 5, model_type = "dem", directed = TRUE), "must be named")
  
  arg_bad_name <- list(data = list(a = mat))
  attr(arg_bad_name, "process") <- "0-1"
  expect_error(redeem:::InitRedeemTerm.dyadic_cov(arg_bad_name, 5, model_type = "dem", directed = TRUE), "must be numeric")
  
  arg_bad_elem <- list(data = list("0" = 1:2))
  attr(arg_bad_elem, "process") <- "0-1"
  expect_error(redeem:::InitRedeemTerm.dyadic_cov(arg_bad_elem, 5, model_type = "dem", directed = TRUE), "must be a matrix or a scalar")
  
  arg_bad_dim <- list(data = list("0" = matrix(0, 4, 4)))
  attr(arg_bad_dim, "process") <- "0-1"
  expect_error(redeem:::InitRedeemTerm.dyadic_cov(arg_bad_dim, 5, model_type = "dem", directed = TRUE), "must be a 5 x 5 matrix")
  
  arg_no_zero <- list(data = list("10" = mat))
  attr(arg_no_zero, "process") <- "0-1"
  expect_error(redeem:::InitRedeemTerm.dyadic_cov(arg_no_zero, 5, model_type = "dem", directed = TRUE), "include at least one measurement time <= 0")
})

test_that("monadic_cov coverage", {
  arg_base <- list()
  attr(arg_base, "process") <- "0-1"

  # Missing data/fun
  expect_error(redeem:::InitRedeemTerm.monadic_cov(arg_base, 5, model_type = "dem", directed = TRUE), "requires a 'data' argument")
  
  arg_no_fun <- list(data = 1:5)
  attr(arg_no_fun, "process") <- "0-1"
  expect_error(redeem:::InitRedeemTerm.monadic_cov(arg_no_fun, 5, model_type = "dem", directed = TRUE), "requires a 'fun' argument")
  
  # Static monadic
  arg_static <- list(data = 1:5, fun = `+`)
  attr(arg_static, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm.monadic_cov(arg_static, 5, model_type = "dem", directed = FALSE), "list")
  expect_type(redeem:::InitRedeemTerm.monadic_cov(arg_static, 5, model_type = "dem", directed = TRUE), "list")
  
  arg_short <- list(data = 1:4, fun = `+`)
  attr(arg_short, "process") <- "0-1"
  expect_error(redeem:::InitRedeemTerm.monadic_cov(arg_short, 5, model_type = "dem", directed = FALSE), "does not match n_nodes")

  
  # Time-varying monadic
  l_mon <- list("0" = 1:5, "10" = 1:5)
  arg_tv <- list(data = l_mon, fun = `+`)
  attr(arg_tv, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm.monadic_cov(arg_tv, 5, model_type = "dem", directed = FALSE), "list")
  expect_type(redeem:::InitRedeemTerm.monadic_cov(arg_tv, 5, model_type = "dem", directed = TRUE), "list")
  
  arg_tv_short <- list(data = list("0" = 1:4), fun = `+`)
  attr(arg_tv_short, "process") <- "0-1"
  expect_error(redeem:::InitRedeemTerm.monadic_cov(arg_tv_short, 5, model_type = "dem", directed = FALSE), "does not match n_nodes")
})

test_that("baseline coverage", {
  arg_base <- list()
  attr(arg_base, "process") <- "0-1"

  # Missing changepoints
  expect_error(redeem:::InitRedeemTerm.baseline(arg_base, 5, model_type = "dem", directed = TRUE), "requires a 'changepoints' argument")
  
  arg_bad_cp <- list(changepoints = "string")
  attr(arg_bad_cp, "process") <- "0-1"
  expect_error(redeem:::InitRedeemTerm.baseline(arg_bad_cp, 5, model_type = "dem", directed = TRUE), "must be numeric")
  
  # Valid
  arg_ok <- list(changepoints = c(10, 20))
  attr(arg_ok, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm.baseline(arg_ok, 5, model_type = "dem", directed = TRUE), "list")
  
  arg_ok2 <- list(changepoint = c(10, 20))
  attr(arg_ok2, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm.baseline(arg_ok2, 5, model_type = "dem", directed = TRUE), "list")
  
  # Labels
  arg_lab <- list(changepoints = c(20, 10), labels = c("B", "A"))
  attr(arg_lab, "process") <- "0-1"
  expect_type(redeem:::InitRedeemTerm.baseline(arg_lab, 5, model_type = "dem", directed = TRUE), "list")
  
  arg_lab_err <- list(changepoints = c(10, 20), labels = c("A"))
  attr(arg_lab_err, "process") <- "0-1"
  expect_error(redeem:::InitRedeemTerm.baseline(arg_lab_err, 5, model_type = "dem", directed = TRUE), "same length")
})

