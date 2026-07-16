# Core Estimation Logic for REM and DEM Transitions

This internal helper function encapsulates the estimation routines for a
single relational or durational event transition. It is used by both
\`rem()\` and \`dem()\`.

## Usage

``` r
estimate_transition(
  data,
  formula_original,
  formula_new,
  indicators,
  n_nodes,
  estimate_method,
  it_max,
  tol,
  accelerated,
  subsample,
  verbose,
  estimate_degree,
  directed,
  semiparametric = FALSE,
  labels_changepoints = NULL,
  time_changepoints = NULL,
  coef_init = NULL,
  model_type = "dem",
  process = "0-1",
  return_data = TRUE,
  save_hist = TRUE,
  use_glm = FALSE,
  legacy = FALSE,
  inf_unidentifiable = TRUE,
  events = NULL
)
```
