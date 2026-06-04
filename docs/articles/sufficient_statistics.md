# Mathematical Definitions of Sufficient Statistics

## Overview

The **redeem** package models the intensities of interaction formation
and dissolution using a log-linear formulation. The intensity (rate) for
a dyad \\(i,j)\\ at time \\t\\ is:

\\\lambda\_{i,j}(t) = \exp(s\_{i,j}(\mathscr{H}\_t)^\top \beta +
\alpha_i + \alpha_j + f(t, \gamma))\\

where: - \\s\_{i,j}(\mathscr{H}\_t)\\ is a vector of **sufficient
statistics** derived from the history of interactions \\\mathscr{H}\_t\\
up to time \\t\\. - \\\alpha_i, \alpha_j\\ are actor popularity
parameters. - \\f(t, \gamma) = \sum\_{q=1}^Q \gamma_q
\mathbb{I}(c\_{q-1} \le t \< c_q)\\ is the baseline step-function
representing temporal variation over change points \\0 = c_0 \< c_1 \<
\dots \< c_Q\\ (with \\\gamma_1 = 0\\).

This vignette provides precise mathematical definitions for all
sufficient network statistics implemented in the package, as well as a
guide for adding custom statistics.

The complete list of available formula terms is documented in the
`redeem_terms` reference manual page. For model estimation, see the help
pages for the models via `dem` and `rem`.

------------------------------------------------------------------------

## Statistic Transformations

Each sufficient statistic \\s\_{i,j}(\mathscr{H}\_t)\\ is defined as a
transformed count \\f(N(t))\\, where \\N(t)\\ is a raw network
statistic. The package supports five standard transformations
\\f(\cdot)\\:

| Transformation | Mathematical Definition       |
|:---------------|:------------------------------|
| **`identity`** | \\f(s) = s\\                  |
| **`log`**      | \\f(s) = \log(1 + s)\\        |
| **`recip`**    | \\f(s) = 1/(1+s)\\            |
| **`bin`**      | \\f(s) = \mathbb{I}(s \> 0)\\ |
| **`sig`**      | \\f(s) = \frac{s}{s + K}\\    |

To achieve maximum performance during likelihood estimation and
simulation, the transformations are computed **incrementally**
on-the-fly inside the C++ state matrix. When a raw statistic increases
by a value \\v\\:

- **Identity**: \\f(s + v) = f(s) + v\\.
- **Log**: \\f(s + v) = f(s) + \log(1 + v e^{-f(s)})\\.
- **Reciprocal**: \\f(s + v) = \frac{1}{\frac{1}{f(s)} + v}\\.
- **Binary**: \\f(s + v) = \mathbb{I}(s + v \> 0)\\.
- **Sigmoid**: \\f(s + v) = \frac{K \cdot f(s) + v (1 - f(s))}{K + v
  (1 - f(s))}\\.

------------------------------------------------------------------------

## Endogenous Network Statistics

These statistics capture structural dependencies within the network
evolution.

### 1. Intercept (`Intercept` / `intercept`)

A constant term representing the baseline log-intensity.
\\s\_{i,j}(\mathscr{H}\_t) = 1\\

### 2. Inertia (`inertia` / `number_interaction`)

Counts how many times the dyad \\(i,j)\\ has initiated an interaction in
the past. \\s\_{i,j}(\mathscr{H}\_t) = f(N\_{i,j}(t))\\ where
\\N\_{i,j}(t) = \sum\_{k: t_k \< t} \mathbb{I}(i_k = i, j_k = j)\\ (or
the windowed version \\N\_{i,j}^w(t) = \sum\_{k: t-w \< t_k \< t}
\mathbb{I}(i_k = i, j_k = j)\\).

### 3. Reciprocity (`reciprocity`)

Models the tendency to reciprocate past interactions (directed only).
\\s\_{i,j}(\mathscr{H}\_t) = f(N\_{j,i}(t))\\

### 4. Duration (`duration` / `current_interaction`)

Measures the dependency on the time since the current interaction
started (DEM only). \\s\_{i,j}(\mathscr{H}\_t) = \begin{cases} f(t -
t\_{\text{start, } i,j}) & \text{if dyad } (i,j) \text{ is interacting
at } t \\ 0 & \text{otherwise} \end{cases}\\ where \\t\_{\text{start, }
i,j}\\ is the timestamp of the last formation event.

### 5. Participation Shifts (P-shifts)

P-shifts capture the sequential dependencies between consecutive events
(REM only). Let the preceding event in the stream be \\A \to B\\. For a
candidate event \\i \to j\\ at time \\t\\:

- **`psABBA`** (Reciprocation): \\s\_{i,j}(t) = \mathbb{I}(i = B, j =
  A)\\
- **`psABBY`** (Receiver turn-continuing): \\s\_{i,j}(t) = \mathbb{I}(i
  = B, j \ne A, j \ne B)\\
- **`psABAY`** (Sender turn-continuing): \\s\_{i,j}(t) = \mathbb{I}(i =
  A, j \ne A, j \ne B)\\
- **`psABXA`** (Usurpation to sender): \\s\_{i,j}(t) = \mathbb{I}(i \ne
  A, i \ne B, j = A)\\
- **`psABXB`** (Usurpation to receiver): \\s\_{i,j}(t) = \mathbb{I}(i
  \ne A, i \ne B, j = B)\\
- **`psABXY`** (Completely new dyad): \\s\_{i,j}(t) = \mathbb{I}(i \ne
  A, i \ne B, j \ne A, j \ne B)\\

------------------------------------------------------------------------

## Triadic Closure and Shared Partners

Triadic statistics capture structural closure. They can be calculated
over active edges (designated as **`current_`** in DEM) or historical
event existence (designated as **`general_`** in REM/DEM). Let
\\\mathcal{A}\_t\\ represent the set of interacting dyads (for
`current_`) or previously interacted dyads (for `general_`) at time
\\t\\.

### Common Partners (`general_common_partners` / `current_common_partners`)

Counts the number of third-party actors \\k\\ sharing a connection of a
specified type with both \\i\\ and \\j\\:

- **`OSP`** (Outgoing Shared Partner): \\s\_{i,j}(t) = f(\|\\k : (i,k)
  \in \mathcal{A}\_t \land (j,k) \in \mathcal{A}\_t\\\|)\\
- **`ISP`** (Incoming Shared Partner): \\s\_{i,j}(t) = f(\|\\k : (k,i)
  \in \mathcal{A}\_t \land (k,j) \in \mathcal{A}\_t\\\|)\\
- **`OTP`** (Outgoing Two-Path): \\s\_{i,j}(t) = f(\|\\k : (i,k) \in
  \mathcal{A}\_t \land (k,j) \in \mathcal{A}\_t\\\|)\\
- **`ITP`** (Incoming Two-Path): \\s\_{i,j}(t) = f(\|\\k : (k,i) \in
  \mathcal{A}\_t \land (j,k) \in \mathcal{A}\_t\\\|)\\

### Triangles (`general_triangle` / `current_triangle`)

Similar to common partners, but only non-zero if the focal dyad itself
is active (directed only):

- **`OSP`** Triangle: \\s\_{i,j}(t) = \mathbb{I}((i,j) \in
  \mathcal{A}\_t) \times f(\|\\k : (i,k) \in \mathcal{A}\_t \land (j,k)
  \in \mathcal{A}\_t\\\|)\\
- **`ISP`** Triangle: \\s\_{i,j}(t) = \mathbb{I}((i,j) \in
  \mathcal{A}\_t) \times f(\|\\k : (k,i) \in \mathcal{A}\_t \land (k,j)
  \in \mathcal{A}\_t\\\|)\\
- **`OTP`** Triangle: \\s\_{i,j}(t) = \mathbb{I}((i,j) \in
  \mathcal{A}\_t) \times f(\|\\k : (i,k) \in \mathcal{A}\_t \land (k,j)
  \in \mathcal{A}\_t\\\|)\\
- **`ITP`** Triangle: \\s\_{i,j}(t) = \mathbb{I}((i,j) \in
  \mathcal{A}\_t) \times f(\|\\k : (k,i) \in \mathcal{A}\_t \land (j,k)
  \in \mathcal{A}\_t\\\|)\\

------------------------------------------------------------------------

## Degree and Centrality Statistics

These statistics capture actor-level activity or popularity. Let
\\\mathcal{D}\_t\\ be the network state at time \\t\\, and \\d\_{i,
\text{out}}(t)\\, \\d\_{i, \text{in}}(t)\\ represent the out-degree and
in-degree of \\i\\ in \\\mathcal{D}\_t\\. Let \\c\_{i, \text{out}}(t)\\,
\\c\_{i, \text{in}}(t)\\ be the total out-events and in-events involving
\\i\\.

### Degree Statistics (`general_degree_out_sender`, etc.)

- **Out-Degree Sender**: \\s\_{i,j}(t) = f(d\_{i, \text{out}}(t))\\
- **Out-Degree Receiver**: \\s\_{i,j}(t) = f(d\_{j, \text{out}}(t))\\
- **In-Degree Sender**: \\s\_{i,j}(t) = f(d\_{i, \text{in}}(t))\\
- **In-Degree Receiver**: \\s\_{i,j}(t) = f(d\_{j, \text{in}}(t))\\
- **Degree Sum** (undirected only): \\s\_{i,j}(t) = f(d\_{i}(t) +
  d\_{j}(t))\\
- **Degree Absolute Difference** (undirected only): \\s\_{i,j}(t) =
  f(\\d\_{i}(t) - d\_{j}(t)\\)\\

### Count Statistics (`general_count_out_sender`, etc.)

Count statistics are identical to degree statistics but use total
interaction counts (\\c\\) rather than binary degrees (\\d\\): -
**Out-Count Sender**: \\s\_{i,j}(t) = f(c\_{i, \text{out}}(t))\\ -
**Out-Count Receiver**: \\s\_{i,j}(t) = f(c\_{j, \text{out}}(t))\\ -
**In-Count Sender**: \\s\_{i,j}(t) = f(c\_{i, \text{in}}(t))\\ -
**In-Count Receiver**: \\s\_{i,j}(t) = f(c\_{j, \text{in}}(t))\\ -
**Count Sum** (undirected only): \\s\_{i,j}(t) = f(c\_{i}(t) +
c\_{j}(t))\\ - **Count Absolute Difference** (undirected only):
\\s\_{i,j}(t) = f(\\c\_{i}(t) - c\_{j}(t)\\)\\

All of these degree and count statistics can be conveniently specified
in model formulas using the
[`degree()`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
(or
[`degrees()`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md))
and
[`count()`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
helper functions:

- **Out-Degree / Out-Count Sender**: `degree(type = "out_sender")` or
  `count(type = "out_sender")`
- **Out-Degree / Out-Count Receiver**: `degree(type = "out_receiver")`
  or `count(type = "out_receiver")`
- **In-Degree / In-Count Sender**: `degree(type = "in_sender")` or
  `count(type = "in_sender")`
- **In-Degree / In-Count Receiver**: `degree(type = "in_receiver")` or
  `count(type = "in_receiver")`
- **Degree / Count Sum**: `degree(type = "sum")` or
  `count(type = "sum")`
- **Degree / Count Absolute Difference**: `degree(type = "absdiff")` or
  `count(type = "absdiff")`

------------------------------------------------------------------------

## Exogenous Statistics

### 1. Dyadic Covariate (`dyadic_cov`)

\\s\_{i,j}(\mathscr{H}\_t) = f(X\_{i,j}(t))\\ where \\X(t)\\ is an \\N
\times N\\ matrix.

### 2. Monadic Covariate (`monadic_cov`)

\\s\_{i,j}(\mathscr{H}\_t) = g(x_i(t), x_j(t))\\ where \\x\\ is a
monadic vector and \\g(\cdot)\\ is a user-defined function.

------------------------------------------------------------------------

## Developer Guide: Adding Custom Sufficient Statistics

To add a new sufficient statistic (e.g., `my_stat`) to the package,
follow this two-step process:

### Step 1: Implement the Update Logic in C++

Define a C++ function in the `src/` directory (e.g., in a new file or in
`src/sufficient_statistics.cpp`) with the signature `ValidateFunction`.
This function computes the change to the statistic state matrix when an
event occurs:

``` cpp
#include "redeem/sufficient_statistics.h"
#include "redeem/extension_api.hpp"

arma::uvec stat_my_stat(
    Data_DEM &object, 
    arma::mat &data, 
    unsigned int &from, 
    unsigned int &to, 
    unsigned int &number_event, 
    unsigned int col_number, 
    std::string transformation, 
    unsigned int type
) {
  if (from == 0 || to == 0) return arma::uvec();
  
  // 1. Calculate which dyad indices are affected and the change value (val)
  double val = (type == 1) ? 1.0 : -1.0; 
  arma::uvec affected_indices = { object.current_stats.find_from_to(from, to) };
  
  // 2. Apply the update using the helper (handles transformations internally)
  apply_update(object, affected_indices, col_number, val, transformation, data);
  
  // 3. Return the indices of modified dyads
  return affected_indices;
}

// 4. Register the term in the global C++ Registry
TERM_REGISTER("my_stat", stat_my_stat);
```

### Step 2: Write the R Initializer Function

In `R/init_terms.R`, create an initialization function named
`InitRedeemTerm.my_stat`. This function validates user-provided
arguments and sets up initial values:

``` r

#' @keywords internal
InitRedeemTerm.my_stat <- function(arglist, n_nodes, model_type, directed, ...) {
  # Validate arguments against expected types and models
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  
  # Define initial statistic values at time 0 (optional)
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  
  # Return the configuration list
  list(
    base_name = "my_stat",
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream
  )
}
```

Once these steps are completed, run
[`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
to update namespaces and documentation, and build the package. The new
term is now ready for use in
[`dem()`](https://corneliusfritz.github.io/redeem/reference/dem.md) and
[`rem()`](https://corneliusfritz.github.io/redeem/reference/rem.md) by
specifying it directly in the formula (e.g.,
`~ my_stat(transformation = "log")`).
