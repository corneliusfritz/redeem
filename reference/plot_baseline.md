# Plot the Estimated Baseline Intensity

Draws a step-function plot of the estimated piecewise-constant baseline
intensity against time. The function dispatches to class-specific
methods for
[`dem`](https://corneliusfritz.github.io/redeem/reference/dem.md),
[`rem`](https://corneliusfritz.github.io/redeem/reference/rem.md), and
`redeem_result` objects.

## Usage

``` r
plot_baseline(x, ...)
```

## Arguments

- x:

  A [`dem`](https://corneliusfritz.github.io/redeem/reference/dem.md),
  [`rem`](https://corneliusfritz.github.io/redeem/reference/rem.md), or
  `redeem_result` object produced by
  [`dem`](https://corneliusfritz.github.io/redeem/reference/dem.md) or
  [`rem`](https://corneliusfritz.github.io/redeem/reference/rem.md).

- ...:

  Additional arguments passed to
  [`graphics::plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

The original object `x` is returned invisibly. Called primarily for its
side effect of producing a plot.
