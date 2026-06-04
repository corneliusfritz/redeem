# Predict the baseline intensity trend at one or more future time points

Decomposes the estimated piecewise-constant log-baseline (`est_time`)
into a smooth trend component (via LOESS) and a seasonal/residual
component, following the same approach used in the application plot
script. The fitted trend is then *extrapolated* to `target_times` using
`predict.loess()` so that the baseline used for out-of-sample scoring
reflects the long-run level of activity rather than any arbitrary fixed
value.

## Usage

``` r
predict_baseline_trend(model, target_times, loess_span = 0.75)
```

## Arguments

- model:

  A `redeem_result` object with non-null `est_time` and
  `time_changepoints` fields.

- target_times:

  Numeric vector: the times at which to predict the trend (typically the
  unique timestamps of the test events).

- loess_span:

  Numeric; the span argument passed to
  [`stats::loess`](https://rdrr.io/r/stats/loess.html). Larger values
  give a smoother (more conservative) trend extrapolation. Defaults to
  0.75.

## Value

A numeric vector of predicted log-baseline (trend component) values, one
per element of `target_times`. Falls back to `mean(est_time)` for each
time point if there are fewer than 3 observations or if the LOESS fit
fails.

## Details

The decomposition mirrors the plot code in the application:

1.  Build a data frame of `(time, est_time)` using the changepoints
    stored in `model$time_changepoints`. The first interval \[0,
    changepoint_1) is given time = 0; each subsequent interval gets the
    corresponding changepoint value.

2.  Fit LOESS on the log-scale `est_time` values.

3.  Predict at each `target_times`; predictions are clamped to the range
    of the observed `est_time` to avoid wild extrapolation.
