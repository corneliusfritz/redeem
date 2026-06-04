# Relational Event Model (REM) Estimation

Estimates a Relational Event Model (REM) for network data, focusing on
the incidence of discrete events between pairs of actors. See
[`dem`](https://corneliusfritz.github.io/redeem/reference/dem.md) for
the full Durational Event Model, which extends the REM to handle
interactions with non-negligible duration.

## Usage

``` r
rem(
  events,
  training_start = 0,
  exogenous_end = NULL,
  formula = NULL,
  n_nodes = NULL,
  directed = FALSE,
  coef = NULL,
  semiparametric = FALSE,
  control = control.redeem()
)
```

## Arguments

- events:

  A matrix of events with columns `time`, `from`, `to`, and optionally
  `type` (1 for start, 3 for exogenous changes).

- training_start:

  Numeric; the time point at which to start the estimation. Defaults to
  0.

- exogenous_end:

  Numeric; optional end time for exogenous baseline changes. Defaults to
  NULL.

- formula:

  A one-sided [`formula`](https://rdrr.io/r/stats/formula.html)
  specifying the sufficient statistics to include in the intensity
  function. The right-hand side must be composed of terms from
  [`redeem_terms`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md).
  For example: `~ inertia() + reciprocity() + degree`. An intercept
  (`~ 1`) is the minimal specification. Defaults to NULL.

- n_nodes:

  Integer; the total number of actors in the network. If `NULL`
  (default), it is automatically identified based on the actors in the
  `events` set.

- directed:

  Logical; whether the interaction events are directed. Defaults to
  FALSE.

- coef:

  Numeric vector; initial coefficients for the model. If provided, this
  must be a concatenated vector of:

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

- control:

  A list of control parameters from
  [`control.redeem`](https://corneliusfritz.github.io/redeem/reference/control.redeem.md).
  Defaults to
  [`control.redeem()`](https://corneliusfritz.github.io/redeem/reference/control.redeem.md).

## Value

An object of class
[`rem_object`](https://corneliusfritz.github.io/redeem/reference/rem_object.md)
containing model estimates and log-likelihoods. See
[`rem_object`](https://corneliusfritz.github.io/redeem/reference/rem_object.md)
for details on the components of the returned object and S3 methods.

## Details

The REM can be viewed as the incidence sub-model of the full
[`dem`](https://corneliusfritz.github.io/redeem/reference/dem.md),
corresponding to the formation process \\\lambda^{0\rightarrow 1}\\. It
uses a counting process approach to estimate the influence of various
covariates on the timing and occurrence of events, assuming that events
are instantaneous points in time.

## Model Formulation

The Relational Event Model characterizes the instantaneous rate at which
actor pair \\(i,j)\\ initiates an event. Under the log-linear
specification, the event intensity at time \\t\\ is:
\$\$\lambda\_{i,j}(t \mid \mathscr{H}\_t,\\ \theta) =
\exp\\\bigl(s\_{i,j}(\mathscr{H}\_t)^\top \alpha + \beta_i + \beta_j +
f(t, \gamma)\bigr)\$\$ where:

- \\s\_{i,j}(\mathscr{H}\_t)\\ is a vector of sufficient statistics
  computed from the event history \\\mathscr{H}\_t\\; see
  [`redeem_terms`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  for available terms.

- \\\alpha\\ is the vector of covariate effects.

- \\\beta_i\\ and \\\beta_j\\ are optional actor-specific baselines
  (sender and receiver sociality), included via the bare symbol `degree`
  in the formula.

- \\f(t, \gamma)\\ is an optional piecewise-constant temporal baseline,
  included via `baseline(changepoints)` in the formula.

## Semiparametric Baseline

When `semiparametric = TRUE`, the temporal baseline rate of event
occurrence is left completely unspecified, and the model parameters are
estimated via the Cox partial likelihood using the `survival` package.
In this path:

- Each observed event time is treated as a failure time, and all
  non-occurring dyads at that time constitute the risk set.

- The exact waiting times between events are conditioned away, meaning
  that inference is based solely on the sequence of events and the
  relative dyadic intensities.

- This approach is equivalent to the *ordered* (or *conditional*) REM
  likelihood introduced by Butts (2008). It is highly robust to temporal
  fluctuations and baseline misspecification since no piecewise baseline
  or changepoints need to be specified.

- **Limitations**: This path does **not** support the specialized
  scalable estimation of sender/receiver popularity effects (`degree`)
  or piecewise-constant temporal baselines.

## References

Fritz, C., Rastelli, R., Fop, M., & Caimo, A. (2026). Scalable
Durational Event Models: Application to Physical and Digital
Interactions. arXiv:2504.00049.

Butts, C. T. (2008). A Relational Event Framework for Social Action.
Sociological Methodology, 38(1), 155-200.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate some relational event data
n <- 20
events <- matrix(c(
  1.2, 1, 5,
  3.1, 2, 8,
  4.5, 1, 3
), ncol = 3, byrow = TRUE)
colnames(events) <- c("time", "from", "to")

# Estimate a simple REM
fit <- rem(
  events = events,
  n_nodes = n,
  formula = ~1,
  control = control.redeem(it_max = 50)
)
summary(fit)
} # }
```
