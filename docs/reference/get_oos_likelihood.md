# Out-of-sample Log-Likelihood (Proper Scoring Rule)

This function computes the out-of-sample log-likelihood (a strictly
proper scoring rule) for each test event under a fitted REM or DEM.

## Usage

``` r
get_oos_likelihood(
  object,
  verbose = FALSE,
  edgelist_test,
  edgelist_train = NULL,
  baseline_method = c("last", "trend", "mean", "beginning"),
  loess_span = 0.75
)
```

## Arguments

- object:

  A `redeem` object (either
  [`rem`](https://corneliusfritz.github.io/redeem/reference/rem.md) or
  [`dem`](https://corneliusfritz.github.io/redeem/reference/dem.md)).

- verbose:

  Logical; if \`TRUE\`, prints verbose output. Defaults to FALSE.

- edgelist_test:

  A matrix or data frame of test events (timing, from, to, type).

- edgelist_train:

  A matrix or data frame of train events (timing, from, to, type).
  Defaults to \`NULL\`, in which case it retrieves the training events
  from the \`object\` or the preprocessed data.

- baseline_method:

  Character; how to compute the fixed log-baseline intensity used for
  out-of-sample scoring. One of: \`"last"\` (uses the last estimated
  baseline value), \`"trend"\` (extrapolates a LOESS trend), \`"mean"\`,
  or \`"beginning"\`. Defaults to \`"last"\`.

- loess_span:

  Numeric; LOESS span (0, 1\] passed to
  [`predict_baseline_trend`](https://corneliusfritz.github.io/redeem/reference/predict_baseline_trend.md)
  when `baseline_method = "trend"`. Defaults to 0.75.

## Value

A numeric vector of log-likelihoods for each test event.

## See also

[`rem_object`](https://corneliusfritz.github.io/redeem/reference/rem_object.md)
and
[`dem_object`](https://corneliusfritz.github.io/redeem/reference/dem_object.md)
for details on prediction methods.
