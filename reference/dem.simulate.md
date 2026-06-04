# Simulate events based on specified formulas and coefficients

Simulate events based on specified formulas and coefficients

## Usage

``` r
dem.simulate(
  events = matrix(0, 0, 4),
  formula_0_1 = NULL,
  formula_1_0 = NULL,
  coef_0_1 = numeric(0),
  coef_1_0 = numeric(0),
  coef_degree_0_1 = 0,
  coef_degree_1_0 = 0,
  n_events = 0,
  time = 0,
  max_events = 4e+05,
  n_nodes,
  verbose = FALSE,
  baseline_0_1 = NULL,
  baseline_1_0 = NULL,
  simultaneous_interactions = TRUE,
  seed = 123,
  directed = FALSE
)
```

## Arguments

- events:

  A matrix representing the initial events with columns `time`, `from`,
  `to`, and `type` (1 for start, 0 for end, 3 for exogenous changes).
  Defaults to an empty 4-column matrix.

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

- coef_0_1:

  Numeric vector; coefficients for the formation process (\\0
  \rightarrow 1\\). Defaults to an empty numeric vector.

- coef_1_0:

  Numeric vector; coefficients for the dissolution process (\\1
  \rightarrow 0\\). Defaults to an empty numeric vector.

- coef_degree_0_1:

  Numeric; degree coefficient for the formation process (\\0 \rightarrow
  1\\). Defaults to 0.

- coef_degree_1_0:

  Numeric; degree coefficient for the dissolution process (\\1
  \rightarrow 0\\). Defaults to 0.

- n_events:

  Integer; number of events to simulate. Defaults to 0.

- time:

  Numeric; simulation time limit. Defaults to 0.

- max_events:

  Integer; maximum number of total events. Defaults to 400000.

- n_nodes:

  Integer; the total number of actors in the network.

- verbose:

  Logical; if `TRUE`, prints progress information. Defaults to FALSE.

- baseline_0_1:

  Numeric vector; baseline for the 0 to 1 transition. If the formula for
  this process contains an
  [`Intercept`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  or a
  [`degree`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  term, then `baseline_0_1` should be a numeric vector with length equal
  to the number of changepoints, representing the shifts in the baseline
  for each interval after the first. If the formula contains neither,
  then `baseline_0_1` must have length equal to the number of
  changepoints + 1. Defaults to NULL.

- baseline_1_0:

  Numeric vector; baseline for the 1 to 0 transition. Similar to
  `baseline_0_1`, its length depends on whether the 1 to 0 formula
  contains an intercept or degree term. Defaults to NULL.

- simultaneous_interactions:

  Logical; whether to allow simultaneous interactions (i.e. multiple
  active events for the same actor or dyad at the same time). Defaults
  to TRUE.

- seed:

  Integer; random seed for simulation. Defaults to 123.

- directed:

  Logical; whether the interaction events are directed. Defaults to
  FALSE.

## Value

A matrix of simulated events.

## Note

Multi-stream event models are currently not supported in simulation.

## Examples

``` r
# Simulate events from a DEM model structure
n <- 10
f_0_1 <- ~ 1 + inertia(transformation = "identity")
f_1_0 <- ~ 1

# Simulating events
evs <- dem.simulate(
  formula_0_1 = f_0_1,
  formula_1_0 = f_1_0,
  n_nodes = n,
  time = 2.0,
  coef_0_1 = c(1.0, 0.5),
  coef_1_0 = c(-0.5),
  seed = 42,
  max_events = 100
)
head(evs)
#>            [,1] [,2] [,3] [,4]
#> [1,] 0.02256350    5   10    1
#> [2,] 0.05274257    2    7    1
#> [3,] 0.06051135    5    9    1
#> [4,] 0.07809648    3    5    1
#> [5,] 0.07831044    1    7    1
#> [6,] 0.09898228    2    6    1
```
