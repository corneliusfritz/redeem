# Control Parameters for REDEEM Models

Unified control object to manage estimation parameters for
[`rem`](https://corneliusfritz.github.io/redeem/reference/rem.md) and
[`dem`](https://corneliusfritz.github.io/redeem/reference/dem.md)
functions.

## Usage

``` r
control.redeem(
  it_max = 100,
  tol = 1e-10,
  accelerated = FALSE,
  verbose = FALSE,
  weighting = TRUE,
  subsample = 1,
  build_time = NULL,
  use_glm = FALSE,
  return_data = FALSE,
  save_hist = TRUE,
  estimate = "Blockwise",
  legacy = FALSE,
  check_matrix = FALSE,
  inf_unidentifiable = TRUE,
  mpl = FALSE
)
```

## Arguments

- it_max:

  Integer; maximum number of iterations for the algorithm. Defaults to
  100.

- tol:

  Numeric; convergence tolerance. Defaults to 1e-10.

- accelerated:

  Logical; if `TRUE`, uses SQUAREM acceleration for MM updates. Defaults
  to FALSE.

- verbose:

  Logical; if `TRUE`, prints progress information. Defaults to FALSE.

- weighting:

  Logical; whether to use weighting to group identical observations.
  Defaults to TRUE.

- subsample:

  Numeric; proportion of data to subsample for internal GLM checks.
  Defaults to 1.

- build_time:

  Numeric; time at which to start building the estimation dataset.
  Events before this time are used to compute statistics but not
  included as observations. Defaults to NULL, in which case all events
  are included.

- use_glm:

  Logical; if `TRUE`, uses standard GLM for updating core coefficients.
  This is often slower but can yield more robust updates. Defaults to
  FALSE.

- return_data:

  Logical; whether to return preprocessed data frames in the result.
  Defaults to FALSE.

- save_hist:

  Logical; whether to save the iteration history of coefficients.
  Defaults to TRUE.

- estimate:

  Character; estimation method for
  [`dem`](https://corneliusfritz.github.io/redeem/reference/dem.md) and
  [`rem`](https://corneliusfritz.github.io/redeem/reference/rem.md)
  ("Blockwise", "NR", or "GD"). Defaults to "Blockwise".

- legacy:

  Logical; if `TRUE`, uses a single `glm.fit` call instead of the
  iterative loop. Defaults to FALSE.

- check_matrix:

  Logical; whether to apply
  [`check_matrix`](https://corneliusfritz.github.io/redeem/reference/check_matrix.md)
  to the event data before estimation. If `TRUE`, repairs missing events
  (e.g., adding start events for interactions that only have end
  events). Defaults to FALSE.

- inf_unidentifiable:

  Logical; whether to set unidentifiable coefficients (e.g., actors with
  0 event counts, globally invariant/collinear covariates) to `-Inf`.
  Defaults to TRUE.

## Value

A list of class `"redeem_control"` containing the specified parameters.
