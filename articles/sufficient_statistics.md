# Mathematical Definitions and Formula Specifications of Sufficient Statistics

## Overview

The **redeem** package provides a flexible, log-linear framework for
modeling interaction event sequences in relational networks. Models are
fitted using
[`rem()`](https://corneliusfritz.github.io/redeem/reference/rem.md) for
**Relational Event Models (REM)** and
[`dem()`](https://corneliusfritz.github.io/redeem/reference/dem.md) for
**Durational Event Models (DEM)**.

In both model classes, event intensities are specified on the right-hand
side of R formulas passed to
[`rem()`](https://corneliusfritz.github.io/redeem/reference/rem.md) or
[`dem()`](https://corneliusfritz.github.io/redeem/reference/dem.md) via
structural statistics \\s\_{i,j}(\mathscr{H}\_t)\\ calculated from event
history \\\mathscr{H}\_t\\:

\\\lambda\_{i,j}(t) = \exp\left(s\_{i,j}(\mathscr{H}\_t)^\top \beta +
\alpha_i + \alpha_j + f(t, \gamma)\right)\\

where:

- \\s\_{i,j}(\mathscr{H}\_t)\\ is a vector of **sufficient statistics**
  derived from history \\\mathscr{H}\_t\\ up to time \\t\\.
- \\\alpha_i, \alpha_j\\ are actor popularity / activity fixed effects
  (specified via `degree` or `degrees`).
- \\f(t, \gamma) = \sum\_{q=1}^Q \gamma_q \mathbb{I}(c\_{q-1} \le t \<
  c_q)\\ is a baseline step-function over change points \\0 = c_0 \< c_1
  \< \dots \< c_Q\\ (specified via
  [`baseline()`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)).

Terms are specified as function calls on the right-hand side of model
formulas (e.g., `formula = ~ intercept() + inertia() + reciprocity()`).
This vignette details the mathematical definitions of all supported
sufficient statistics alongside examples of how they are specified in
formula objects.

For details on estimation procedures and options, see the help pages for
[`rem()`](https://corneliusfritz.github.io/redeem/reference/rem.md),
[`dem()`](https://corneliusfritz.github.io/redeem/reference/dem.md), and
[`control.redeem()`](https://corneliusfritz.github.io/redeem/reference/control.redeem.md).

------------------------------------------------------------------------

## Statistic Transformations

Each sufficient statistic \\s\_{i,j}(\mathscr{H}\_t)\\ can be modified
by applying a transformation function \\f(\cdot)\\ to a raw network
statistic \\x\\. The package supports five standard transformations:

| Transformation | Mathematical Definition | R Formula Syntax |
|:---|:---|:---|
| **`identity`** (default) | \\f(x) = x\\ | `transformation = "identity"` |
| **`log`** | \\f(x) = \log(1 + x)\\ | `transformation = "log"` |
| **`recip`** | \\f(x) = \frac{1}{1 + x}\\ | `transformation = "recip"` |
| **`bin`** | \\f(x) = \mathbb{I}(x \> 0)\\ | `transformation = "bin"` |
| **`sig`** | \\f(x) = \frac{x}{x + K}\\ | `transformation = "sig", K = 1` |

#### Formula Examples

Transformations can be passed to any statistic via the `transformation`
argument (and optional `K` parameter):

``` r

# Default identity transformation
formula_1 <- ~ inertia(transformation = "identity")
# Log-transformed count statistic
formula_2 <- ~ inertia(transformation = "log")
# Binary indicator (has dyad interacted at least once?)
formula_3 <- ~ reciprocity(transformation = "bin")
# Sigmoid-like saturation with parameter K = 5
formula_4 <- ~ general_common_partners(type = "OSP", transformation = "sig", K = 5)
```

------------------------------------------------------------------------

## Baseline and Nuisance Terms

These terms model baseline rates, step-function temporal variations, and
actor-level fixed effects.

#### 1. Intercept (`intercept` / `Intercept`)

Represents a constant baseline log-intensity across all dyads.

\\s\_{i,j}(\mathscr{H}\_t) = 1\\

``` r

~ intercept()
```

#### 2. Time-Varying Baseline (`baseline`)

Captures temporal variation by fitting a piece-wise constant baseline
across user-specified change points \\0 = c_0 \< c_1 \< \dots \< c_Q\\.

\\s\_{i,j}(\mathscr{H}\_t) = \sum\_{q=1}^Q \mathbb{I}(c\_{q-1} \le t \<
c_q)\\

``` r

# Piece-wise constant baseline changing at times 10, 25, and 50
~ baseline(changepoints = c(10, 25, 50))
# Baseline with custom labels
~ baseline(changepoints = c(100, 200), labels = c("Phase1", "Phase2", "Phase3"))
```

#### 3. Degree Fixed Effects (`degree` / `degrees`)

Adds node-specific sender and receiver baseline parameters (\\\alpha_i\\
and \\\gamma_j\\) estimated efficiently via Minorization-Maximization
(MM) steps. For directed events, fixed sender and receiver effects are
included: \\\alpha_i + \gamma_j\\ To make the model identifiable, we set
the first finite sender parameter to zero. For undirected events, fixed
popularity effects are included: \\\alpha_i + \alpha_j\\ In this case,
no constraints are needed for identifiability. However, no additional
intercept can be included. For isolated actors, the corresponding fixed
effect is set to \\-\infty\\.

``` r

~ degrees + reciprocity()
```

------------------------------------------------------------------------

## Endogenous Dyadic Terms

Endogenous dyadic terms capture direct past interaction tendencies
between unit pairs \\(i,j)\\.

#### 1. Inertia (`inertia` / `number_interaction`)

Measures the tendency for dyads with past interactions to interact
again. Let \\N\_{i,j}(t)\\ be the total events from \\i\\ to \\j\\
before time \\t\\, and \\N\_{i,j}^w(t) = \sum\_{k: t-w \< t_k \< t}
\mathbb{I}(i_k = i, j_k = j)\\ be the count within time window \\w\\.

\\s\_{i,j}(\mathscr{H}\_t) = f(N\_{i,j}(t)) \quad \text{or} \quad
f(N\_{i,j}^w(t))\\

``` r

# Cumulative interaction count (identity transformation)
~ inertia()
# Log-transformed inertia evaluated over a 30-day window
~ inertia(transformation = "log", window = 30)
# Inertia calculated from a secondary event stream
~ inertia(event_stream = external_events)
# Alternative name for inertia
~ number_interaction(transformation = "log")
```

#### 2. Reciprocity (`reciprocity`)

Measures the tendency of \\i \to j\\ interactions to be triggered by
past \\j \to i\\ interactions (directed networks only).

\\s\_{i,j}(\mathscr{H}\_t) = f(N\_{j,i}(t)) \quad \text{or} \quad
f(N\_{j,i}^w(t))\\

``` r

# Cumulative reciprocity
~ reciprocity()
# Binary reciprocity indicator over a 14-unit window
~ reciprocity(transformation = "bin", window = 14)
```

#### 3. Duration (`duration` / `current_interaction`)

Measures dependency on the time elapsed since the current interaction
started (DEM dissolution process only).

\\s\_{i,j}(\mathscr{H}\_t) = \begin{cases} f(t - t\_{\text{start}, i,j})
& \text{if dyad } (i,j) \text{ is active at } t \\ 0 & \text{otherwise}
\end{cases}\\

``` r

# Elapsed duration (identity transformation)
~ duration()
# Log-transformed duration for DEM dissolution formulas
~ duration(transformation = "log")
# Alternative name
~ current_interaction(transformation = "log")
```

#### 4. Participation Shifts (`ps` / `psABBA`, `psABBY`, `psABAY`, `psABXA`, `psABXB`, `psABXY`)

Participation shifts capture immediate sequential dependencies between
consecutive events (REM only). Let the previous event be \\A \to B\\.
For candidate event \\C \to D\\ at time \\t\\:

| Term Function | Shorthand Mode | Math Definition | Description |
|:---|:---|:---|:---|
| [`psABBA()`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md) | `mode = "ABBA"` | \\\mathbb{I}(C = B, D = A)\\ | Direct reciprocation (\\B \to A\\) |
| [`psABBY()`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md) | `mode = "ABBY"` | \\\mathbb{I}(C = B, D \ne A, D \ne B)\\ | Receiver turn-continuing (\\B \to Y\\) |
| [`psABAY()`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md) | `mode = "ABAY"` | \\\mathbb{I}(C = A, D \ne A, D \ne B)\\ | Sender turn-continuing (\\A \to Y\\) |
| [`psABXA()`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md) | `mode = "ABXA"` | \\\mathbb{I}(C \ne A, C \ne B, D = A)\\ | Usurpation to sender (\\X \to A\\) |
| [`psABXB()`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md) | `mode = "ABXB"` | \\\mathbb{I}(C \ne A, C \ne B, D = B)\\ | Usurpation to receiver (\\X \to B\\) |
| [`psABXY()`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md) | `mode = "ABXY"` | \\\mathbb{I}(C \ne A, C \ne B, D \ne A, D \ne B)\\ | Completely new dyad (\\X \to Y\\) |

``` r

~ psABBA() + psABBY() + psABAY()
# Convenience terms
~ ps(mode = "ABBA") + ps(mode = "ABAY")
```

------------------------------------------------------------------------

## Triadic Closure and Shared Partners

Triadic terms capture structural closure via shared third-party nodes
\\k\\. Statistics can be based on historical event existence
(`general_`, REM/DEM) or currently active edges (`current_`, DEM only).

#### 1. Common Partners (`general_common_partners` / `current_common_partners` / `common_partner`)

Counts the number of third-party actors \\k\\ sharing directed paths of
a specified type with \\i\\ and \\j\\:

- **Outgoing Shared Partners (`OSP`)**: Both \\i\\ and \\j\\ interact
  with \\k\\ (\\i \to k\\ and \\j \to k\\). \\s\_{i,j}(t) =
  f\left(\left\|\\k : N\_{i,k}(t)\>0 \land
  N\_{j,k}(t)\>0\\\right\|\right)\\
- **Incoming Shared Partners (`ISP`)**: Both \\i\\ and \\j\\ receive
  interactions from \\k\\ (\\k \to i\\ and \\k \to j\\). \\s\_{i,j}(t) =
  f\left(\left\|\\k : N\_{k,i}(t)\>0 \land
  N\_{k,j}(t)\>0\\\right\|\right)\\
- **Outgoing Two-Path (`OTP`)**: Path from \\i\\ through \\k\\ to \\j\\
  (\\i \to k\\ and \\k \to j\\). \\s\_{i,j}(t) = f\left(\left\|\\k :
  N\_{i,k}(t)\>0 \land N\_{k,j}(t)\>0\\\right\|\right)\\
- **Incoming Two-Path (`ITP`)**: Path from \\j\\ through \\k\\ to \\i\\
  (\\k \to i\\ and \\j \to k\\). \\s\_{i,j}(t) = f\left(\left\|\\k :
  N\_{k,i}(t)\>0 \land N\_{j,k}(t)\>0\\\right\|\right)\\

``` r

# Full term names
~ general_common_partners(type = "OTP", transformation = "log")
~ current_common_partners(type = "OSP")
# Convenience terms
~ common_partner(history = "general", type = "OTP", transformation = "log")
~ common_partner(history = "current", type = "OSP")
```

#### 2. Triangles (`general_triangle` / `current_triangle` / `triangle`)

Similar to common partners, but conditional on the focal dyad \\(i,j)\\
having an active or historical edge (directed networks only):

\\s\_{i,j}(t) = \mathbb{I}((i,j) \in \mathcal{A}\_t) \times
f\left(\|CP\_{i,j}^{\text{type}}(t)\|\right)\\

``` r

# Historical OTP triangle (log transformed)
~ general_triangle(type = "OTP", transformation = "log")
# Active OSP triangle in DEM
~ current_triangle(type = "OSP")
# Convenience terms
~ triangle(history = "general", type = "OTP")
~ triangle(history = "current", type = "OSP")
```

------------------------------------------------------------------------

## Degree and Centrality Statistics

Degree and event-count statistics capture actor activity and popularity.
Let \\d\_{i,\text{out}}(t)\\ and \\d\_{i,\text{in}}(t)\\ be node \\i\\’s
out- and in-degrees at time \\t\\, and \\c\_{i,\text{out}}(t)\\,
\\c\_{i,\text{in}}(t)\\ be node \\i\\’s total sent and received event
counts.

#### 1. Binary Degree Statistics (`general_degree_*`, `current_degree_*`, `degree`)

Captures structural centrality based on unique interaction partners
(binary degree).

- **Sender Out-Degree**: \\s\_{i,j}(t) = f(d\_{i,\text{out}}(t))\\
- **Receiver Out-Degree**: \\s\_{i,j}(t) = f(d\_{j,\text{out}}(t))\\
- **Sender In-Degree**: \\s\_{i,j}(t) = f(d\_{i,\text{in}}(t))\\
- **Receiver In-Degree**: \\s\_{i,j}(t) = f(d\_{j,\text{in}}(t))\\
- **Degree Sum** (undirected networks): \\s\_{i,j}(t) = f(d\_{i}(t) +
  d\_{j}(t))\\
- **Degree Absolute Difference** (undirected networks): \\s\_{i,j}(t) =
  f(\|d\_{i}(t) - d\_{j}(t)\|)\\

``` r

~ general_degree_out_sender() + general_degree_in_receiver()
~ current_degree_sum(transformation = "log")
# Convenience terms
~ degree(type = "out_sender", history = "general")
~ degree(type = "in_receiver", history = "current")
~ degree(type = "sum", history = "general", transformation = "log")
~ degree(type = "absdiff", history = "general")
```

#### 2. Weighted Event Count Statistics (`general_count_*`, `current_count_*`, `count`)

Identical to degree statistics, but computed using total interaction
event counts (\\c\\) rather than binary degrees (\\d\\).

- **Sender Out-Count**: \\s\_{i,j}(t) = f(c\_{i,\text{out}}(t))\\
- **Receiver In-Count**: \\s\_{i,j}(t) = f(c\_{j,\text{in}}(t))\\
- **Count Sum** (undirected networks): \\s\_{i,j}(t) = f(c\_{i}(t) +
  c\_{j}(t))\\
- **Count Absolute Difference** (undirected networks): \\s\_{i,j}(t) =
  f(\|c\_{i}(t) - c\_{j}(t)\|)\\

``` r

~ general_count_out_sender(transformation = "log")
~ current_count_sum()
# Convenience terms
~ count(type = "out_sender", history = "general", transformation = "log")
~ degree(type = "out_sender", history = "general", count = TRUE)
```

------------------------------------------------------------------------

## Exogenous Statistics

Exogenous statistics incorporate external nodal or dyadic covariates
into model formulas.

#### 1. Dyadic Covariate (`dyadic_cov`)

Incorporates an external \\N \times N\\ dyadic covariate matrix \\X(t)\\
(constant or time-varying).

\\s\_{i,j}(\mathscr{H}\_t) = f(X\_{i,j}(t))\\

``` r

# Constant dyadic matrix (must of size N x N)
~ dyadic_cov(data = dist_matrix)
# Time-varying dyadic matrix list across change points
cov_list <- list("0" = matrix1, "50" = matrix2)
~ dyadic_cov(data = cov_list, change_points = c(0, 50))
```

#### 2. Monadic Covariate (`monadic_cov`)

Transforms external node-level attributes \\x(t)\\ into a dyadic matrix
using a user-specified function \\g(x_i, x_j)\\.

\\s\_{i,j}(\mathscr{H}\_t) = g(x_i(t), x_j(t))\\

``` r

# Absolute difference in age (homophily)
~ monadic_cov(data = age_vector, fun = function(u, v) abs(u - v))
# Sender main effect
~ monadic_cov(data = status_vector, fun = function(u, v) u)
# Receiver main effect 
~ monadic_cov(data = status_vector, fun = function(u, v) v)
# Product interaction
~ monadic_cov(data = status_vector, fun = function(u, v) u * v)
```
