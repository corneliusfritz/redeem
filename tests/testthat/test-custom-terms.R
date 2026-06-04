test_that("Dynamic C++ Custom Term Registration and Estimation works perfectly", {
  skip_on_cran()

  # Ensure the package is loaded so C-callables are visible
  library(redeem)
  library(Rcpp)

  # Find include directories manually to bypass compiler-specific depends plugins
  rcpp_include <- system.file("include", package = "Rcpp")
  rcpp_armadillo_include <- system.file("include", package = "RcppArmadillo")
  # Try installed include directory first, then fallback to relative development paths
  inst_include <- system.file("include", package = "redeem")
  if (inst_include == "" || !dir.exists(inst_include)) {
    paths <- c(
      file.path(getwd(), "../../inst/include"),
      file.path(getwd(), "../inst/include"),
      file.path(getwd(), "inst/include")
    )
    for (p in paths) {
      if (dir.exists(p)) {
        inst_include <- normalizePath(p)
        break
      }
    }
  }

  old_flags <- Sys.getenv("PKG_CXXFLAGS")
  old_libs <- Sys.getenv("PKG_LIBS")

  # Add all includes manually and link only to standard BLAS/LAPACK
  Sys.setenv(PKG_CXXFLAGS = paste0(
    "-I", shQuote(inst_include),
    " -I", shQuote(rcpp_include),
    " -I", shQuote(rcpp_armadillo_include)
  ))
  Sys.setenv(PKG_LIBS = "-lRlapack -lRblas")

  on.exit({
    Sys.setenv(PKG_CXXFLAGS = old_flags)
    Sys.setenv(PKG_LIBS = old_libs)
  }, add = TRUE)

  # Define and compile a custom C++ statistic using Rcpp::sourceCpp
  cpp_code <- '
  #include <RcppArmadillo.h>
  #include <redeem/extension_api.hpp>
  #include <redeem/sufficient_statistics.h>

  arma::uvec custom_intercept(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
    if (from == 0 || to == 0) return arma::uvec();
    arma::uvec indices = arma::regspace<arma::uvec>(0, object.current_stats.data.n_rows - 1);
    apply_update(object, indices, col_number, 1.0, transformation, data);
    return indices;
  }

  // [[Rcpp::export]]
  void register_custom_term() {
    typedef void (*reg_fn_t)(const char*, void*);
    reg_fn_t reg = (reg_fn_t)R_GetCCallable("redeem", "redeem_register_term_C");
    if (!reg) {
      Rcpp::stop("Could not retrieve C-callable: redeem_register_term_C");
    }
    reg("custom_intercept", (void*)&custom_intercept);
  }
  '

  # Compile the code in the current session
  Rcpp::sourceCpp(code = cpp_code, verbose = FALSE)

  # Register our custom intercept term!
  register_custom_term()

  # Define R-side term initializer matching the exact signature
  InitRedeemTerm.custom_intercept <- function(arglist, n_nodes, model_type, directed, ...) {
    return(list(
      base_name = "custom_intercept",
      eval_at_zero = matrix(1, n_nodes, n_nodes)
    ))
  }
  assign("InitRedeemTerm.custom_intercept", InitRedeemTerm.custom_intercept, envir = .GlobalEnv)
  on.exit(rm("InitRedeemTerm.custom_intercept", envir = .GlobalEnv), add = TRUE)

  # Prepare simple dummy event matrix
  ed <- matrix(c(
    1, 1, 2,
    2, 2, 3,
    3, 1, 3,
    4, 3, 2,
    5, 2, 1
  ), ncol = 3, byrow = TRUE)
  colnames(ed) <- c("time", "from", "to")

  # Fit a REM model using the custom intercept term!
  fit_custom <- rem(
    events = ed,
    formula = ~ custom_intercept(),
    n_nodes = 3,
    directed = TRUE
  )

  # Fit a REM model using the standard intercept term to verify equality
  fit_standard <- rem(
    events = ed,
    formula = ~ Intercept(),
    n_nodes = 3,
    directed = TRUE
  )

  # Check that coefficients and likelihood are identical!
  expect_equal(fit_custom$coefficients, fit_standard$coefficients, tolerance = 1e-6)
  expect_equal(fit_custom$log_likelihood, fit_standard$log_likelihood, tolerance = 1e-6)
})
