# redeem Model Terms

The help pages of
[`rem`](https://corneliusfritz.github.io/redeem/reference/rem.md) and
[`dem`](https://corneliusfritz.github.io/redeem/reference/dem.md)
describe the model formulation and estimation details. This page
documents all statistics available to be used in the model formulas,
characterizing the intensities of event formation and dissolution.

In the `redeem` framework, models like DEM (fitted via
[`dem`](https://corneliusfritz.github.io/redeem/reference/dem.md)) and
REM (fitted via
[`rem`](https://corneliusfritz.github.io/redeem/reference/rem.md)) are
specified using R formulas. The right-hand side of these formulas
defines the structural statistics and covariates, where each term must
be specified separately as an explicit function call (e.g.,
`~ inertia() + reciprocity(window = 10)`).

All terms support an optional `transformation` argument \\f\\. The
available transformations are:

- `"identity"` (default): \\f(x) = x\\

- `"log"`: \\f(x) = \log(x + 1)\\

- `"recip"`: \\f(x) = 1/(x+1)\\

- `"bin"`: \\f(x) = I(x \> 0)\\

- `"sig"`: sigmoid-like saturation, \\f(x) = x/(x + K)\\

Throughout, \\N\_{i,j}(t)\\ denotes the cumulative number of events from
\\i\\ to \\j\\ up to (but not including) time \\t\\; \\N\_{i,j}^w(t)\\
is the windowed analogue on \\(t-w,\\t)\\; \\d_i^{\mathrm{out}}(t) =
\|\\l: N\_{i,l}(t)\>0\\\|\\ is the historical out-degree of \\i\\; and
\\c_i^{\mathrm{out}}(t) = \sum_l N\_{i,l}(t)\\ is the total event count
sent by \\i\\. The superscript \\\mathrm{act}\\ indicates that the
quantity is computed on the currently active DEM network.

The implemented terms are grouped into five categories:

1.  **Baseline and Nuisance Terms**: Intercept, time-varying baseline,
    and degree fixed effects.

2.  **Endogenous Dyadic Terms**: Inertia, reciprocity, interaction
    duration, and participation shifts.

3.  **Triadic Closure and Shared Partners**: Common partners and
    triangle statistics.

4.  **Degree and Centrality Statistics**: Actor degree and event count
    statistics.

5.  **Exogenous Covariates**: Dyadic and monadic covariate terms.

## Arguments

- K:

  Numeric; the evaluation point or scaling/saturation factor for the
  sufficient statistic (default is 1).

- transformation:

  Character; specifies the transformation to apply to the statistic. One
  of:

  - `"identity"` (default): \\f(x) = x\\

  - `"log"`: \\f(x) = \log(x+1)\\

  - `"recip"`: \\f(x) = 1/(x+1)\\

  - `"bin"`: \\f(x) = I(x\>0)\\

  - `"sig"`: sigmoid-like saturation, \\f(x) = x/(x+K)\\

- event_stream:

  Optional matrix or data frame; an alternative event stream to use for
  calculating the statistic. If `NULL` (default), the modeled stream is
  used.

- window:

  Numeric; time window for calculating the statistic (default `Inf`,
  i.e., use full history).

- type:

  Character; the specific variation of the statistic or triangle type
  (e.g., `"OSP"`, `"ISP"`, `"OTP"`, `"ITP"`, `"out_sender"`, `"sum"`).

- mode:

  Character; the participation shift mode (e.g., `"ABBA"`, `"ABBY"`).

- data:

  For `dyadic_cov`, a numeric matrix of dimensions \\N \times N\\, a
  scalar applied globally, or a named list of matrices for time-varying
  covariates. For `monadic_cov`, a numeric vector of length \\N\\ or a
  named list of vectors for time-varying covariates.

- fun:

  A function taking two arguments `fun(v_i, v_j)` to generate dyadic
  values.

- change_points:

  Optional numeric vector; time points for time-varying covariates if
  `data` is a list.

- changepoints:

  Numeric vector; time points where the baseline intensity is allowed to
  change.

- labels:

  Character vector; optional labels for the resulting time intervals.

- history:

  Character; `"general"` for cumulative history or `"current"` for
  currently active events.

- count:

  Logical; if `TRUE`, uses count-based (weighted) versions of degree
  statistics (default `FALSE`).

- ...:

  Arguments passed to the underlying initialization function.

## Value

A `redeem_term` object (a list containing structural information about
the statistic) to be used inside model formulas.

## Multi-Stream Event Covariates

Most endogenous terms support covariates calculated from multiple event
streams. By providing an `event_stream` argument to a term (e.g.,
`inertia(event_stream = other_events)`), users can model one event
process while accounting for the history of another. The package
automatically handles the splintering and union of these timelines.

## Baseline and Nuisance Terms

- `Intercept()`: Intercept: Constant log-intensity baseline.
  \\s\_{i,j}(t) = 1\\. Also available as `intercept()`.

- `baseline(changepoints, labels)`: Baseline: Stepwise constant
  log-baseline with user-specified change points \\c_1 \< c_2 \< \ldots
  \< c_K\\. \\s\_{i,j}(t) = \sum\_{k=1}^{K} I(t \in \[c_k, c\_{k+1}))\\.
  Coefficients are treated as nuisance parameters.

- `degree` / `degrees`: Degree Fixed Effects: Node-specific sender and
  receiver baselines \\\alpha_i\\ and \\\gamma_j\\ estimated via
  Minorization-Maximization (MM). Contribution to linear predictor:
  \\\alpha_i + \gamma_j\\. Treated as nuisance parameters.

## Endogenous Dyadic Terms

- `inertia(transformation, K, event_stream, window)`: Inertia:
  Cumulative count of past events from \\i\\ to \\j\\. \\s\_{i,j}(t) =
  f(N\_{i,j}(t))\\; windowed: \\f(N\_{i,j}^w(t))\\.

- `reciprocity(transformation, K, event_stream, window)`: Reciprocity:
  Cumulative count of past events from \\j\\ to \\i\\. \\s\_{i,j}(t) =
  f(N\_{j,i}(t))\\ **(Directed only)**.

- `current_interaction(transformation, K, event_stream)`: Duration: Time
  elapsed since the currently active event \\(i,j)\\ started.
  \\s\_{i,j}(t) = f(t - t\_{\mathrm{start},i,j})\\ **(DEM only)**.
  Alias: `duration()`.

- **Participation Shifts** (for two consecutive events \\(A \to B) \to
  (C \to D)\\, each statistic is 1 if the specified pattern holds, 0
  otherwise; **REM only, Directed only**):

  - `psABBA(event_stream)`: PS-ABBA: Receiver responds to sender.
    \\s\_{C,D}(t) = I(C = B,\\ D = A)\\.

  - `psABBY(event_stream)`: PS-ABBY: Receiver initiates to a new target.
    \\s\_{C,D}(t) = I(C = B,\\ D \ne A)\\.

  - `psABAY(event_stream)`: PS-ABAY: Sender initiates to a new target.
    \\s\_{C,D}(t) = I(C = A,\\ D \ne B)\\.

  - `psABXA(event_stream)`: PS-ABXA: Outsider targets original sender.
    \\s\_{C,D}(t) = I(C \ne A, C \ne B,\\ D = A)\\.

  - `psABXB(event_stream)`: PS-ABXB: Outsider targets original receiver.
    \\s\_{C,D}(t) = I(C \ne A, C \ne B,\\ D = B)\\.

  - `psABXY(event_stream)`: PS-ABXY: Entirely new dyad. \\s\_{C,D}(t) =
    I(C \ne A, C \ne B,\\ D \ne A, D \ne B)\\.

  - `ps(mode, event_stream)`: PS Shorthand: Dispatches to one of the six
    participation shift statistics above based on `mode` (one of
    `"ABBA"`, `"ABBY"`, `"ABAY"`, `"ABXA"`, `"ABXB"`, `"ABXY"`).

## Triadic Closure and Shared Partners

- `general_common_partners(` `transformation, K, type,`
  `event_stream, window)`: Historical Common Partners: Number of nodes
  \\k\\ sharing a historical directed path of the specified type with
  both \\i\\ and \\j\\. \\s\_{i,j}(t) =
  f(\|CP\_{i,j}^{\mathrm{type}}(t)\|)\\.

  - `"OSP"` (Outgoing Shared Partner): \\N\_{i,k}(t)\>0\\ and
    \\N\_{j,k}(t)\>0\\.

  - `"ISP"` (Incoming Shared Partner): \\N\_{k,i}(t)\>0\\ and
    \\N\_{k,j}(t)\>0\\.

  - `"OTP"` (Outgoing Two-Path): \\N\_{i,k}(t)\>0\\ and
    \\N\_{k,j}(t)\>0\\.

  - `"ITP"` (Incoming Two-Path): \\N\_{k,i}(t)\>0\\ and
    \\N\_{j,k}(t)\>0\\.

  Aliases: `general_common_partner()`, `general_common_partner_OSP()`,
  `general_common_partner_ISP()`, `general_common_partner_OTP()`,
  `general_common_partner_ITP()`.

- `current_common_partners(` `transformation, K, type,` `event_stream)`:
  Active Common Partners: As `general_common_partners` but restricted to
  currently active edges. \\s\_{i,j}(t) =
  f(\|CP\_{i,j}^{\mathrm{type,act}}(t)\|)\\ **(DEM only)**. Aliases:
  `current_common_partner()`, `current_common_partner_OSP()`,
  `current_common_partner_ISP()`, `current_common_partner_OTP()`,
  `current_common_partner_ITP()`.

- `general_triangle(transformation, K, type, event_stream, window)`:
  Historical Triangles: Number of closed triads around \\(i,j)\\ in the
  historical event network of the specified type. \\s\_{i,j}(t) =
  f(\|\Delta\_{i,j}^{\mathrm{type}}(t)\|)\\ **(Directed only)**.

- `current_triangle(transformation, K, type, event_stream)`: Active
  Triangles: As `general_triangle` but restricted to currently active
  edges. **(DEM only, Directed only)**.

- `common_partner(history, type, ...)`: Common Partner Shorthand:
  Dispatches to `general_common_partners()` (`history="general"`) or
  `current_common_partners()` (`history="current"`).

- `triangle(history, type, ...)`: Triangle Shorthand: Dispatches to
  `general_triangle()` (`history="general"`) or `current_triangle()`
  (`history="current"`).

## Degree and Centrality Statistics

- `general_degree_out_sender(`
  `transformation, K, event_stream, window)`: Sender Out-Degree:
  Historical out-degree of sender \\i\\. \\s\_{i,j}(t) =
  f(d_i^{\mathrm{out}}(t))\\ **(Directed only)**.

- `general_degree_out_receiver(`
  `transformation, K, event_stream, window)`: Receiver Out-Degree:
  Historical out-degree of receiver \\j\\. \\s\_{i,j}(t) =
  f(d_j^{\mathrm{out}}(t))\\ **(Directed only)**.

- `general_degree_in_sender(`
  `transformation, K, event_stream, window)`: Sender In-Degree:
  Historical in-degree of sender \\i\\. \\s\_{i,j}(t) =
  f(d_i^{\mathrm{in}}(t))\\ **(Directed only)**.

- `general_degree_in_receiver(`
  `transformation, K, event_stream, window)`: Receiver In-Degree:
  Historical in-degree of receiver \\j\\. \\s\_{i,j}(t) =
  f(d_j^{\mathrm{in}}(t))\\ **(Directed only)**.

- `general_degree_sum(transformation, K, event_stream, window)`: Degree
  Sum: Sum of historical degrees of both endpoints. \\s\_{i,j}(t) =
  f(d_i(t) + d_j(t))\\ **(Undirected only)**.

- `general_degree_absdiff(` `transformation, K, event_stream, window)`:
  Degree Absolute Difference: Absolute difference in historical degrees.
  \\s\_{i,j}(t) = f(\|d_i(t) - d_j(t)\|)\\ **(Undirected only)**.

- `general_count_out_sender(`
  `transformation, K, event_stream, window)`: Sender Out-Count: Total
  events sent by sender \\i\\. \\s\_{i,j}(t) =
  f(c_i^{\mathrm{out}}(t))\\ **(Directed only)**.

- `general_count_out_receiver(`
  `transformation, K, event_stream, window)`: Receiver Out-Count: Total
  events sent by receiver \\j\\. \\s\_{i,j}(t) =
  f(c_j^{\mathrm{out}}(t))\\ **(Directed only)**.

- `general_count_in_sender(` `transformation, K, event_stream, window)`:
  Sender In-Count: Total events received by sender \\i\\. \\s\_{i,j}(t)
  = f(c_i^{\mathrm{in}}(t))\\ **(Directed only)**.

- `general_count_in_receiver(`
  `transformation, K, event_stream, window)`: Receiver In-Count: Total
  events received by receiver \\j\\. \\s\_{i,j}(t) =
  f(c_j^{\mathrm{in}}(t))\\ **(Directed only)**.

- `general_count_sum(transformation, K, event_stream, window)`: Count
  Sum: Sum of total event counts of both endpoints. \\s\_{i,j}(t) =
  f(c_i(t) + c_j(t))\\ **(Undirected only)**.

- `general_count_absdiff(` `transformation, K, event_stream, window)`:
  Count Absolute Difference: Absolute difference in total event counts.
  \\s\_{i,j}(t) = f(\|c_i(t) - c_j(t)\|)\\ **(Undirected only)**.

- `current_degree_out_sender(transformation, K, event_stream)`: Active
  Sender Out-Degree: Out-degree of \\i\\ in the active DEM network.
  \\s\_{i,j}(t) = f(d_i^{\mathrm{out,act}}(t))\\ **(DEM only, Directed
  only)**.

- `current_degree_out_receiver(` `transformation, K, event_stream)`:
  Active Receiver Out-Degree: Out-degree of \\j\\ in active DEM network.
  \\s\_{i,j}(t) = f(d_j^{\mathrm{out,act}}(t))\\ **(DEM only, Directed
  only)**.

- `current_degree_in_sender(transformation, K, event_stream)`: Active
  Sender In-Degree: In-degree of \\i\\ in the active DEM network.
  \\s\_{i,j}(t) = f(d_i^{\mathrm{in,act}}(t))\\ **(DEM only, Directed
  only)**.

- `current_degree_in_receiver(` `transformation, K, event_stream)`:
  Active Receiver In-Degree: In-degree of \\j\\ in active DEM network.
  \\s\_{i,j}(t) = f(d_j^{\mathrm{in,act}}(t))\\ **(DEM only, Directed
  only)**.

- `current_degree_sum(transformation, K, event_stream)`: Active Degree
  Sum: Sum of active degrees of both endpoints. \\s\_{i,j}(t) =
  f(d_i^{\mathrm{act}}(t) + d_j^{\mathrm{act}}(t))\\ **(DEM only,
  Undirected only)**.

- `current_degree_absdiff(transformation, K, event_stream)`: Active
  Degree Absolute Difference: Absolute difference in active degrees.
  \\s\_{i,j}(t) = f(\|d_i^{\mathrm{act}}(t) - d_j^{\mathrm{act}}(t)\|)\\
  **(DEM only, Undirected only)**.

- `current_count_out_sender(transformation, K, event_stream)`: Active
  Sender Out-Count: Total active events sent by \\i\\. \\s\_{i,j}(t) =
  f(c_i^{\mathrm{out,act}}(t))\\ **(DEM only, Directed only)**.

- `current_count_out_receiver(` `transformation, K, event_stream)`:
  Active Receiver Out-Count: Total active events sent by \\j\\.
  \\s\_{i,j}(t) = f(c_j^{\mathrm{out,act}}(t))\\ **(DEM only, Directed
  only)**.

- `current_count_in_sender(transformation, K, event_stream)`: Active
  Sender In-Count: Total active events received by \\i\\. \\s\_{i,j}(t)
  = f(c_i^{\mathrm{in,act}}(t))\\ **(DEM only, Directed only)**.

- `current_count_in_receiver(` `transformation, K, event_stream)`:
  Active Receiver In-Count: Total active events received by \\j\\.
  \\s\_{i,j}(t) = f(c_j^{\mathrm{in,act}}(t))\\ **(DEM only, Directed
  only)**.

- `current_count_sum(transformation, K, event_stream)`: Active Count
  Sum: Sum of active event counts of both endpoints. \\s\_{i,j}(t) =
  f(c_i^{\mathrm{act}}(t) + c_j^{\mathrm{act}}(t))\\ **(DEM only,
  Undirected only)**.

- `current_count_absdiff(transformation, K, event_stream)`: Active Count
  Absolute Difference: Absolute difference in active event counts.
  \\s\_{i,j}(t) = f(\|c_i^{\mathrm{act}}(t) - c_j^{\mathrm{act}}(t)\|)\\
  **(DEM only, Undirected only)**.

- `degree(`
  `history, type, count, transformation, K, event_stream, window)`:
  Degree Shorthand: Dispatches to the appropriate degree statistic based
  on `history` (`"general"` or `"current"`) and `type` (`"out_sender"`,
  `"out_receiver"`, `"in_sender"`, `"in_receiver"`, `"sum"`,
  `"absdiff"`). Set `count = TRUE` for weighted (count-based) variants.
  Alias: `degrees()`.

- `count(history, type, transformation, K, event_stream, window)`: Count
  Shorthand: Equivalent to `degree(..., count = TRUE)`.

## Exogenous Covariates

- `dyadic_cov(data, change_points)`: Dyadic Covariate: Time-constant or
  time-varying external dyadic covariate matrix \\X\\. \\s\_{i,j}(t) =
  X\_{i,j}(t)\\.

- `monadic_cov(data, fun, change_points)`: Monadic Covariate: External
  monadic covariate vector \\x\\ converted to a dyadic matrix via
  user-supplied function \\g\\. \\s\_{i,j}(t) = g(x_i(t),\\ x_j(t))\\.
