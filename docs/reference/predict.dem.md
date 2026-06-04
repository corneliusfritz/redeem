# Predict Method for Durational Event Models (DEM)

Predicts the intensity, linear predictor, or term contributions for a
fitted DEM model.

## Usage

``` r
# S3 method for class 'dem'
predict(
  object,
  time = NULL,
  type = c("response", "lp", "terms"),
  process = c("both", "formation", "dissolution"),
  ...
)
```

## Arguments

- object:

  A
  [`dem_object`](https://corneliusfritz.github.io/redeem/reference/dem_object.md)
  object.

- time:

  Numeric vector; optional time point(s) at which to predict. Defaults
  to NULL, in which case predictions are returned for all intervals in
  the preprocessed data.

- type:

  Character; the type of prediction: \`"response"\` (instantaneous
  intensity), \`"lp"\` (linear predictor), or \`"terms"\` (contributions
  of each term to the linear predictor). Defaults to \`"response"\`.

- process:

  Character; the transition process to predict: \`"both"\` (returns
  predictions for both formation and dissolution), \`"formation"\` (0
  -\> 1 transition), or \`"dissolution"\` (1 -\> 0 transition). Defaults
  to \`"both"\`.

- ...:

  Additional arguments.

## Value

A data frame containing predictions or term contributions.
