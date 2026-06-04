# Simulate a Relational Event Model (REM)

Simulate a Relational Event Model (REM)

## Usage

``` r
rem.simulate(
  events = matrix(0, 0, 4),
  formula,
  coef = NULL,
  coef_degree = 0,
  n_events = 0,
  time = 0,
  max_events = 4e+05,
  n_nodes,
  verbose = FALSE,
  baseline = NULL,
  seed = 123,
  block = 1,
  directed = FALSE
)
```

## Arguments

- events:

  A matrix representing the initial events with columns `time`, `from`,
  `to`, and optionally `type` (1 for start, 3 for exogenous changes).
  Defaults to an empty 4-column matrix.

- formula:

  A one-sided [`formula`](https://rdrr.io/r/stats/formula.html)
  specifying the sufficient statistics to include in the intensity
  function. The right-hand side must be composed of terms from
  [`redeem_terms`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md).
  For example: `~ inertia() + reciprocity() + degree`. An intercept
  (`~ 1`) is the minimal specification.

- coef:

  Numeric vector; coefficients for the model. Defaults to NULL.

- coef_degree:

  Numeric; degree coefficient. Defaults to 0.

- n_events:

  Integer; number of events to simulate. Defaults to 0.

- time:

  Numeric; simulation time. Defaults to 0.

- max_events:

  Integer; maximum number of events. Defaults to 400000.

- n_nodes:

  Integer; the total number of actors in the network.

- verbose:

  Logical; if `TRUE`, prints progress information. Defaults to FALSE.

- baseline:

  Numeric vector; baseline intensity values for intervals defined by
  changepoints. Defaults to NULL.

- seed:

  Integer; random seed. Defaults to 123.

- block:

  An integer vector of length `n_nodes` indicating the block/group
  assignment for each node, or a single value applied to all nodes.
  Defaults to 1. If multiple blocks are assigned, within-block
  interactions are suppressed (i.e., their dyadic intensities are set to
  0), meaning only events occurring between actors in different blocks
  are simulated.

- directed:

  Logical; whether the interaction events are directed. Defaults to
  FALSE.

## Value

A matrix of simulated events.

## Details

The `block` parameter allows the user to specify a partition of the
nodes into different groups (blocks). When the vector contains more than
one unique block identifier:

- The simulation suppresses all within-block dyad intensities by setting
  them to 0.

- Consequently, only events between nodes belonging to different blocks
  are generated (between-block interactions).

- If all nodes belong to the same block (e.g., if a single value or
  `NULL` is passed), no block-level constraints are applied, and all
  dyads are simulated according to the specified model formula.

## Note

Multi-stream event models are currently not supported in simulation.

## Examples

``` r
# Simulate events from a REM model structure
n <- 10
f1 <- ~ 1 + inertia(transformation = "identity")

# Simulating events
evs <- rem.simulate(
  formula = f1,
  n_nodes = n,
  time = 1.0,
  coef = c(1.0, 0.5),
  seed = 42,
  max_events = 100
)
head(evs)
#>            [,1] [,2] [,3]
#> [1,] 0.02256350    5   10
#> [2,] 0.05180009    2    8
#> [3,] 0.05909043    5    9
#> [4,] 0.07507357    3    5
#> [5,] 0.07526190    1    8
#> [6,] 0.09287949    2    6
```
