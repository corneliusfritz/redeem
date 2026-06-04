# Predict Method for Relational Event Models (REM)

Predicts the intensity, linear predictor, or term contributions for a
fitted REM model.

## Usage

``` r
# S3 method for class 'rem'
predict(object, time = NULL, type = c("response", "lp", "terms"), ...)
```

## Arguments

- object:

  A
  [`rem_object`](https://corneliusfritz.github.io/redeem/reference/rem_object.md)
  object.

- time:

  Numeric vector; optional time point(s) at which to predict. Defaults
  to NULL, in which case predictions are returned for all intervals in
  the preprocessed data.

- type:

  Character; the type of prediction: \`"response"\` (instantaneous
  intensity), \`"lp"\` (linear predictor), or \`"terms"\` (contributions
  of each term to the linear predictor). Defaults to \`"response"\`.

- ...:

  Additional arguments.

## Value

A data frame containing predictions or term contributions.
