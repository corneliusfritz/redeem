# MM Algorithm for Durational Event Models with Time-Varying Effects

Implementation of the scalable block-coordinate ascent algorithm for
DEMs with time-varying baseline intensities. This function performs
iterative Minorization-Maximization (MM) updates for degree and temporal
effects while using Newton-Raphson for core covariates.

## Usage

``` r
estimate_mmt(
  data,
  indicators,
  it_max,
  n_nodes,
  tol = 1e-10,
  accelerated = TRUE,
  time_changepoints = NULL,
  labels_changepoints = NULL,
  subsample = 0.2,
  verbose = FALSE,
  est_degree = NULL,
  est_core = NULL,
  est_time = NULL,
  estimate_degree = TRUE,
  directed = FALSE,
  return_data = TRUE,
  save_hist = TRUE,
  use_glm = FALSE,
  inf_unidentifiable = TRUE
)
```

## Arguments

- data:

  Preprocessed data.table.

- indicators:

  Numeric vector of covariate indices.

- it_max:

  Maximum number of iterations.

- n_nodes:

  Number of nodes.

- tol:

  Convergence tolerance.

- accelerated:

  Logical; use SQUAREM acceleration for degree effects.

- time_changepoints:

  Numeric vector of time changepoints.

- labels_changepoints:

  Character vector of labels for time slices.

- subsample:

  Subsampling rate for GLM backup estimation.

- verbose:

  Logical; print progress.

- est_degree:

  Initial degree coefficients.

- est_core:

  Initial core coefficients.

- est_time:

  Initial time effects.

- estimate_degree:

  Logical; estimate degree effects.

- directed:

  Logical; whether the network is directed.

- return_data:

  Logical; whether to return the preprocessed data in the result.

- save_hist:

  Logical; whether to save the iteration history of coefficients.

- use_glm:

  Logical; whether to use GLM-based core updates as fallback or control.

- inf_unidentifiable:

  Logical; if TRUE, unidentifiable parameters are set to -Inf.

## Details

The algorithm decomposes the log-likelihood and updates blocks of
parameters sequentially. Specifically:

1.  Core effects (\\\beta\\) are updated via Newton-Raphson.

2.  Degree effects (\\\alpha\\) are updated using an MM step that avoids
    explicit Hessian inversion for high-dimensional actor sets.

3.  Temporal effects (\\\gamma\\) are updated via a similar MM step
    across defined time changepoints.
