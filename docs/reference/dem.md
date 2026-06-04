# Durational Event Model (DEM) Estimation

Estimates a Durational Event Model (DEM) for relational event sequences
where interactions have a duration.

## Usage

``` r
dem(
  events,
  training_start = 0,
  exogenous_end = NULL,
  formula_0_1 = NULL,
  formula_1_0 = NULL,
  n_nodes,
  directed = FALSE,
  estimate_0_1 = NULL,
  estimate_1_0 = NULL,
  coef_0_1 = NULL,
  coef_1_0 = NULL,
  semiparametric = FALSE,
  simultaneous_interactions = TRUE,
  control = control.redeem()
)
```

## Arguments

- events:

  A matrix of events with columns `time`, `from`, `to`, and `type` (1
  for start, 0 for end, 3 for exogenous changes).

- training_start:

  Numeric; the time point at which to start the estimation. Defaults to
  0.

- exogenous_end:

  Numeric; the exogenous time point at which the observational period
  ends. Defaults to NULL, which implies that time when the final event
  was observed is taken as the end of the observational period.

- formula_0_1:

  A one-sided [`formula`](https://rdrr.io/r/stats/formula.html)
  specifying the sufficient statistics for the formation process (\\0
  \rightarrow 1\\). The right-hand side must be composed of terms from
  [`redeem_terms`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md).
  For example: `~ inertia() + degree`. An intercept (`~ 1`) is the
  minimal specification. Defaults to NULL.

- formula_1_0:

  A one-sided [`formula`](https://rdrr.io/r/stats/formula.html)
  specifying the sufficient statistics for the dissolution process (\\1
  \rightarrow 0\\). The right-hand side must be composed of terms from
  [`redeem_terms`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md).
  An intercept (`~ 1`) is the minimal specification. Defaults to NULL.

- n_nodes:

  Integer; the total number of actors in the network.

- directed:

  Logical; whether the interaction events are directed. Defaults to
  FALSE.

- estimate_0_1:

  Logical; whether to estimate the formation process. Defaults to NULL,
  in which case it is estimated if `formula_0_1` is provided.

- estimate_1_0:

  Logical; whether to estimate the dissolution process. Defaults to
  NULL, in which case it is estimated if `formula_1_0` is provided.

- coef_0_1:

  Numeric vector; initial coefficients for the formation model. If
  provided, this must be a concatenated vector of:

  1.  Core coefficients: values for sufficient statistics in the
      formula.

  2.  Degree coefficients (if `degree` is in the formula): a vector of
      length `n_nodes` (undirected) or `2 * n_nodes` (directed, sender
      effects first then receiver effects).

  3.  Baseline coefficients (if temporal changepoints are present): a
      vector of length equal to the number of baseline intervals (equal
      to number of changepoints if an intercept/degree is present, or
      changepoints + 1 if neither is present).

  Defaults to NULL, in which case default starting values are
  automatically computed.

- coef_1_0:

  Numeric vector; initial coefficients for the dissolution model. If
  provided, this must be a concatenated vector of:

  1.  Core coefficients: values for sufficient statistics in the
      formula.

  2.  Degree coefficients (if `degree` is in the formula): a vector of
      length `n_nodes` (undirected) or `2 * n_nodes` (directed, sender
      effects first then receiver effects).

  3.  Baseline coefficients (if temporal changepoints are present): a
      vector of length equal to the number of baseline intervals (equal
      to number of changepoints if an intercept/degree is present, or
      changepoints + 1 if neither is present).

  Defaults to NULL, in which case default starting values are
  automatically computed.

- semiparametric:

  Logical; whether to use a semiparametric baseline. Defaults to FALSE.
  See the 'Semiparametric Baseline' section for details.

- simultaneous_interactions:

  Logical; whether to allow simultaneous interactions (i.e. multiple
  active events for the same actor or dyad at the same time). Defaults
  to TRUE.

- control:

  A list of control parameters from
  [`control.redeem`](https://corneliusfritz.github.io/redeem/reference/control.redeem.md).
  Defaults to
  [`control.redeem()`](https://corneliusfritz.github.io/redeem/reference/control.redeem.md).

## Value

An object of class
[`dem_object`](https://corneliusfritz.github.io/redeem/reference/dem_object.md)
containing model estimates, log-likelihoods, and preprocessed data. See
[`dem_object`](https://corneliusfritz.github.io/redeem/reference/dem_object.md)
for details on the components of the returned object and S3 methods.

## Details

The Durational Event Model (DEM) is a general framework for analyzing
durational events, extending standard Relational Event Models (REM) by
decoupling the modeling of event incidence from event duration. It
characterizes the dynamics via two separate continuous-time counting
processes:

- Formation Process (\\0 \rightarrow 1\\):

  Counts the number of times that actor pair \\(i,j)\\ starts an
  interaction up to time \\t\\. The incidence intensity is denoted by
  \\\lambda\_{i,j}^{0\rightarrow 1}(t \| \mathscr{H}\_t)\\.

- Dissolution Process (\\1 \rightarrow 0\\):

  Counts the number of times that actor pair \\(i,j)\\ stops interacting
  up to time \\t\\. The dissolution intensity is denoted by
  \\\lambda\_{i,j}^{1\rightarrow 0}(t \| \mathscr{H}\_t)\\.

Under the assumption that the processes are non-homogeneous Poisson
processes, the intensities are modeled as:
\$\$\lambda\_{i,j}^{0\rightarrow 1}(t \| \mathscr{H}\_t,
\theta^{0\rightarrow 1}) = \exp(s\_{i,j}^{0\rightarrow
1}(\mathscr{H}\_t)^\top \alpha^{0\rightarrow 1} + \beta_i^{0\rightarrow
1} + \beta_j^{0\rightarrow 1} + f(t, \gamma^{0\rightarrow 1}))\$\$
\$\$\lambda\_{i,j}^{1\rightarrow 0}(t \| \mathscr{H}\_t,
\theta^{1\rightarrow 0}) = \exp(s\_{i,j}^{1\rightarrow
0}(\mathscr{H}\_t)^\top \alpha^{1\rightarrow 0} + \beta_i^{1\rightarrow
0} + \beta_j^{1\rightarrow 0} + f(t, \gamma^{1\rightarrow 0}))\$\$
where:

- \\s\_{i,j}(\mathscr{H}\_t)\\ is a vector of dynamic network statistics
  capturing the history of past interactions \\\mathscr{H}\_t\\.

- \\\alpha\\ is a parameter vector determining the covariate effects.

- \\\beta_i\\ and \\\beta_j\\ are actor-specific sociality/popularity
  parameters (degree correction) capturing actor heterogeneity.

- \\f(t, \gamma)\\ is a piecewise-constant step function modeling
  temporal baseline fluctuations across a set of changepoints.

To satisfy the Feller criterion and ensure that the continuous-time
counting process remains non-explosive, count-based network statistics
(such as inertia or common partners) are typically log-transformed on
the \\\log(x + 1)\\ scale.

## Scalable Estimation Algorithm

The likelihood of the model is separable with respect to
\\\theta^{0\rightarrow 1}\\ and \\\theta^{1\rightarrow 0}\\, allowing
independent estimation of the incidence and duration components.
Traditional maximum likelihood estimation via standard Newton-Raphson
requires computing and inverting an \\O(N^2)\\ Hessian matrix, which is
computationally prohibitive for larger networks. To bypass this, the
`redeem` package implements a highly scalable block-coordinate ascent
algorithm that separates parameter updates:

1.  **Step 1**: Update covariate parameters \\\alpha\\ using a standard
    Newton-Raphson update.

2.  **Step 2**: Update high-dimensional actor popularity baselines
    \\\beta\\ using Minorization-Maximization (MM) steps, avoiding
    explicit matrix inversion.

3.  **Step 3**: Update baseline step function parameters \\\gamma\\ via
    a closed-form step.

More information is provided in Fritz et at. (2026).

## Semiparametric Baseline

When `semiparametric = TRUE`, the baseline rates for both the formation
(\\0 \rightarrow 1\\) and dissolution (\\1 \rightarrow 0\\) processes
are left completely unspecified. Both processes are estimated as
separate Cox proportional hazards models using the `survival` package.
In this path:

- For the formation process, each start event (1) is treated as a
  failure, and all inactive dyads at that time constitute the risk set.

- For the dissolution process, each end event (0) is treated as a
  failure, and all currently active interactions constitute the risk
  set.

- The exact waiting times (durations of the active and inactive states)
  are conditioned out, and estimation is based solely on the ordering of
  events and the relative dyadic intensities at each transition time.

- This approach is highly robust to arbitrary temporal fluctuations in
  baseline rates since no piecewise-constant temporal baselines or
  changepoints need to be specified.

- **Limitations**: This path does **not** support the specialized
  scalable estimation of sender/receiver popularity effects (`degree`)
  or piecewise-constant temporal baselines.

## References

Fritz, C., Rastelli, R., Fop, M., & Caimo, A. (2026). Scalable
Durational Event Models: Application to Physical and Digital
Interactions. arXiv:2504.00049.

## Examples

``` r
# Simulate some durational data
n <- 20
events <- matrix(c(
  1.2, 1, 5, 1,
  2.5, 1, 5, 0,
  3.1, 2, 8, 1,
  4.4, 2, 8, 0
), ncol = 4, byrow = TRUE)
colnames(events) <- c("time", "from", "to", "type")

# Estimate a simple DEM
fit <- dem(
  events = events,
  n_nodes = n,
  formula_0_1 = ~1,
  formula_1_0 = ~1,
  control = control.redeem(estimate = "Blockwise")
)
summary(fit)
#> Call:
#> dem(events = events, formula_0_1 = ~1, formula_1_0 = ~1, n_nodes = n, 
#>     control = control.redeem(estimate = "Blockwise"))
#> 
#> Results for Incidence Intensity (0 -> 1): 
#> Fixed Effects:
#>           Estimate Std. Error t value  Pr(>|t|)    
#> Intercept  -6.0324     0.7071 -8.5311 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -14.065 
#> 
#> Results for Duration Intensity (1 -> 0): 
#> Fixed Effects:
#>           Estimate Std. Error t value Pr(>|t|)
#> Intercept -0.26236    0.70711  -0.371   0.7106
#> 
#> Log-likelihood: -2.525 
#> 
#> Combined Model Fit:
#>   AIC: 37.17892 
#>   BIC: 35.95151 
#> 
#> Total estimation time: 0.006180048 secs 
```
