# Standardize Estimation Output for Redeem Models

Standardize Estimation Output for Redeem Models

## Usage

``` r
standardize_redeem_result(
  coefficients,
  coef_hist,
  covariance,
  llh_hist,
  data,
  prediction,
  method = "dem.nr",
  n_nodes = NULL,
  directed = FALSE,
  time_changepoints = NULL,
  labels_changepoints = NULL,
  full_baseline = FALSE,
  return_data = TRUE,
  save_hist = TRUE
)
```

## Arguments

- coefficients:

  Joint coefficient vector.

- coef_hist:

  Matrix of coefficient history.

- covariance:

  Covariance matrix for core effects.

- llh_hist:

  Vector of log-likelihood history.

- data:

  Input data.

- prediction:

  Predicted values.

- method:

  Estimation method name (for class).

- n_nodes:

  Number of nodes.

- directed:

  Logical; are the interaction events directed?

- time_changepoints:

  Numeric vector of time changepoints.

- labels_changepoints:

  Character vector of labels for the changepoints.

- return_data:

  Logical; should the estimation dataset be returned?

- save_hist:

  Logical; should the parameter history be returned?
