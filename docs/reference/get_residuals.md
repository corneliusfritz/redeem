# Get residuals for model diagnostics (Cox-Snell Residuals)

Computes Cox-Snell residuals for a fitted model to diagnose
goodness-of-fit and calibration.

## Usage

``` r
get_residuals(object, get_0_1 = TRUE, get_1_0 = TRUE, raw = FALSE)
```

## Arguments

- object:

  A `redeem` object (either
  [`rem`](https://corneliusfritz.github.io/redeem/reference/rem.md) or
  [`dem`](https://corneliusfritz.github.io/redeem/reference/dem.md)).

- get_0_1:

  Logical; if \`TRUE\`, computes residuals for the formation (0 -\> 1)
  process. Defaults to TRUE.

- get_1_0:

  Logical; if \`TRUE\`, computes residuals for the dissolution (1 -\> 0)
  process. Defaults to TRUE.

- raw:

  Logical; if \`TRUE\`, returns the raw Cox-Snell residuals. Defaults to
  FALSE.

## Value

If \`raw = TRUE\`, a list containing the raw residuals for the selected
process(es). If \`raw = FALSE\`, a list of data frames containing the
Kaplan-Meier coordinates (\`time\`, \`surv\`) and the corresponding
\`theoretical\` standard exponential survival values.

## Details

Cox-Snell residuals are a standard diagnostic tool for continuous-time
survival models and counting processes. Under the true model
specification, the integrated cumulative intensity computed up to the
exact time of an observed event is distributed as a standard exponential
random variable, i.e., \\\Lambda\_{ij}(t_k) \sim Exp(1)\\.

Consequently, if the model is correctly specified:

- The empirical survival function of these residuals should closely
  match the theoretical survival function of a standard exponential
  distribution, \\S(r) = \exp(-r)\\.

- Deviations between the empirical Kaplan-Meier curve of the residuals
  and the theoretical exponential curve signal model misspecification,
  unmodeled dyadic heterogeneity, or non-stationarity.

The function can compute residuals for both the formation/incidence (\\0
\rightarrow 1\\) process and the dissolution/duration (\\1 \rightarrow
0\\) process.

## References

Cox, D. R., & Snell, E. J. (1968). A general definition of residuals.
Journal of the Royal Statistical Society: Series B (Methodological),
30(2), 248-265.
