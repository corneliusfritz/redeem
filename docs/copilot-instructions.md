# Copilot Instructions for redem

## Project Overview

This is an R package implementing **Relational Event Models (REM)** and
**Durational Event Models (DEM)** in a counting process framework. The
package analyzes temporal interaction networks by modeling both the
start (0→1) and end (1→0) of interactions between actors/nodes.

## Architecture

### Core Components

- **R Layer** (`R/estimate_dem.R`): Primary user interface with
  [`dem()`](https://corneliusfritz.github.io/redeem/reference/dem.md)
  and
  [`dem.simulate()`](https://corneliusfritz.github.io/redeem/reference/dem.simulate.md)
  functions. Handles formula preprocessing, data transformation, and
  estimation orchestration.
- **C++ Core** (`src/`): Performance-critical preprocessing and
  computation via Rcpp/RcppArmadillo
  - `DEM` class: Manages durational event models with simultaneous
    interaction support
  - `REM` class: Simplified relational event models (no simultaneous
    interactions)
  - `Current_Stat`: Tracks current state of all actor pairs with O(1)
    lookups using index calculations
  - `Data_DEM`: Container for current/historical statistics and
    interaction tracking
  - `Hist_Events`: Maintains historical interaction records for network
    statistics

### Key Data Flow

1.  User provides event matrix: `[time, from, to, type]` where type=1
    (start) or type=0 (end)
2.  [`formula_preprocess()`](https://corneliusfritz.github.io/redeem/reference/formula_preprocess.md)
    parses R formulas into term lists, transformations, and exogenous
    data
3.  C++ `preprocess()` converts events into long-format survival data
    with sufficient statistics
4.  Estimation via Newton-Raphson or blockwise methods (leveraging
    [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html))

## Critical Patterns

### Index Calculations (Current_Stat)

All pair lookups use **closed-form index formulas** for O(1) access: -
**Directed**: For pair (from, to), index =
`(n_nodes-1)*(from-1) + (to < from ? to-1 : to-2)` - **Undirected**:
Index = `n_nodes*(from-1) - (from-1)*from/2 + (to-from) - 1` (requires
from \< to) - See `find_from_to()`, `find_from()`, `find_to()` in
[current_stat.cpp](https://corneliusfritz.github.io/redeem/src/current_stat.cpp#L8-L23)

### Statistic Transformations

Sufficient statistics support 4 transformations specified as strings: -
`"identity"`: Direct counts (default) - `"log"`: log(exp(x) + val) for
log-scale updates - `"recip"`: 1/(1/x + val) for reciprocal
transformations  
- `"bin"`: Binary indicators (0/1)

Applied via `add_stats()`, `log_add_stats()`, `recip_add_stats()` in
`Current_Stat`. See
[sufficient_statistics.h](https://corneliusfritz.github.io/redeem/src/sufficient_statistics.h#L120-L135)
for usage.

### Formula Syntax

Terms follow pattern: `statistic(data=matrix, transformation=type)` -
Example:
`~ current_interaction + current_common_partners + cov_symm(data=covariates, transformation="log")` -
Parsed by
[`formula_preprocess_single()`](https://corneliusfritz.github.io/redeem/reference/formula_preprocess_single.md)
which extracts data matrices and transformation specs - Two formulas
required: `formula_0_1` (start model) and `formula_1_0` (end model)

### Simultaneous Interactions

When `simultaneous_interactions=FALSE`, actors become unavailable during
interactions: - `not_avail(actor)` sets availability columns (5, 6) to 0
for all pairs involving actor - `now_avail(actor)` restores availability
when interaction ends - Affects risk set in preprocessing

## Build & Development

### Building the Package

``` r
# From R console
devtools::load_all()              # Quick iteration during development
devtools::document()              # Update documentation
devtools::install()               # Full installation

# Or from terminal
R CMD INSTALL --preclean .
```

### Memory Diagnostics

The package uses AddressSanitizer for memory debugging (see
[Makevars](https://corneliusfritz.github.io/redeem/src/Makevars)):

``` make
CXXFLAGS += -fsanitize=address -fno-omit-frame-pointer -g
```

To disable for production builds, comment out the sanitizer flags.

### Rcpp Integration

- All C++ exports defined with `// [[Rcpp::export]]`
- Run
  [`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html)
  after adding/modifying exports
- Auto-generates `RcppExports.{cpp,R}` (do not edit manually)

## Common Tasks

### Adding New Sufficient Statistics

1.  Define function in
    [sufficient_statistics.h](https://corneliusfritz.github.io/redeem/src/sufficient_statistics.h)
    with signature:

    ``` cpp
    arma::uvec stat_name(Data_DEM &object, arma::mat &data, unsigned int &from,
                         unsigned int &to, unsigned int &number_event,
                         unsigned int col_number, std::string transformation, bool type)
    ```

2.  Register in term lookup (implementation-specific)

3.  Return indices of affected pairs; function internally updates
    `object.current_stats.data`

### Debugging Preprocessed Data

Use `return_data=TRUE` in
[`dem()`](https://corneliusfritz.github.io/redeem/reference/dem.md) to
inspect the long-format data:

``` r

result <- dem(events, formula_0_1, formula_1_0, n_nodes=10, return_data=TRUE)
head(result$preprocessed)  # Columns: time_new, time, pair_id, status, event, from, to, from_avail, to_avail, [statistics]
```

## Important Conventions

- **1-based indexing** in R layer; converted to 0-based for internal
  matrix lookups in C++
- **No self-loops**: All validation checks ensure from ≠ to
- Event type: 1 = interaction start (0→1), 0 = interaction end (1→0), 3
  = exogenous baseline change
- Time changepoints: Used for time-varying baseline intensities;
  inserted as type=3 pseudo-events
- Verbose mode (`verbose=TRUE`): Prints progress during preprocessing
  (useful for large datasets)

## Testing & Validation

- No formal test suite currently; validate through simulation:

  ``` r

  sim <- dem.simulate(formula_0_1=~intercept, formula_1_0=~intercept,
                      coef_0_1=c(-2), coef_1_0=c(-1), n_nodes=5, n_events=100)
  fit <- dem(sim, formula_0_1=~intercept, formula_1_0=~intercept, n_nodes=5)
  ```

- Check coefficient recovery and preprocessing correctness
