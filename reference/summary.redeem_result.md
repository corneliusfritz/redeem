# Summary of a `redeem_result` Model Fit

Computes a summary of a fitted `redeem_result` object, collecting the
estimated fixed effects, log-likelihood, and (if present) degree and
temporal baseline effects into a structured list suitable for printing.

## Usage

``` r
# S3 method for class 'redeem_result'
summary(object, ...)
```

## Arguments

- object:

  A `redeem_result` object.

- ...:

  Additional arguments (currently unused).

## Value

An object of class `summary.redeem_result`, which is a list containing:

- `coefficients`: A numeric matrix with one row per fixed-effect
  covariate and columns `Estimate`, `Std. Error`, `t value`, and
  `Pr(>|t|)`.

- `llh`: The log-likelihood of the fitted model (scalar).

- `degree_summary`: A list with summary statistics (`n`,
  `n_unidentifiable`, `mean`, `sd`, `range`) of the estimated degree
  effects, only present when more than 10 degree parameters were
  estimated.

- `degree_effects`: A named numeric vector of estimated degree effects,
  only present when 10 or fewer degree parameters were estimated.

- `time_summary`: A list with summary statistics of the estimated
  temporal baseline effects, only present when more than 10 time
  intervals were used.

- `time_effects`: A named numeric vector of estimated temporal baseline
  effects, only present when 10 or fewer time intervals were used.

- `iter`: Integer; the number of iterations performed by the optimizer
  (`NA` if history was not saved).
