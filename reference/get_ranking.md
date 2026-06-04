# Get ranking for test events (Out-of-Sample Goodness-of-Fit)

Evaluates the out-of-sample predictive performance of a fitted model on
a test event sequence using a ranking-based Goodness-of-Fit (GoF)
procedure.

## Usage

``` r
get_ranking(
  object,
  verbose = FALSE,
  k_max = 1000,
  edgelist_test,
  edgelist_train = NULL,
  ties.method = c("average", "first", "last", "random", "max", "min"),
  return_probabilities = FALSE,
  baseline_method = c("trend", "mean", "last", "beginning"),
  loess_span = 0.75
)
```

## Arguments

- object:

  A `redeem` object (either
  [`rem`](https://corneliusfritz.github.io/redeem/reference/rem.md) or
  [`dem`](https://corneliusfritz.github.io/redeem/reference/dem.md)).

- verbose:

  Logical; if \`TRUE\`, prints verbose output. Defaults to FALSE.

- k_max:

  Maximum number of ranked pairs to return. Defaults to 1000.

- edgelist_test:

  A matrix of test events (timing, from, to, type).

- edgelist_train:

  A matrix of train events (timing, from, to, type). Defaults to NULL.

- ties.method:

  Character; the method to handle ties when ranking event intensities,
  passed directly to [`rank`](https://rdrr.io/r/base/rank.html).
  Defaults to `"average"`. One of:

  `"average"`

  :   Assigns the average of the ranks of all tied elements to each.

  `"first"`

  :   Breaks ties by the order they appear in the data structure.

  `"last"`

  :   Breaks ties by the reverse order of their appearance.

  `"random"`

  :   Breaks ties randomly, ensuring no systematic bias.

  `"max"`

  :   Assigns the maximum of the ranks of the tied elements to all.

  `"min"`

  :   Assigns the minimum of the ranks of the tied elements to all.

- return_probabilities:

  Logical; if TRUE, returns the predicted probabilities/scores instead
  of recall curves. Defaults to FALSE.

- baseline_method:

  Character; how to compute the fixed log-baseline intensity used for
  out-of-sample scoring. Defaults to `"trend"`. One of:

  `"trend"`

  :   Fit a LOESS trend to the estimated piecewise-constant log-baseline
      (`est_time`) over training time and extrapolate it to the start of
      the test period. This mirrors the trend decomposition used in the
      application plot script and typically yields a better forecast
      than a fixed mean.

  `"mean"`

  :   Use the simple mean of `est_time`.

  `"last"`

  :   Use the last estimated log-baseline value (i.e.\\ the value from
      the most recent training interval).

  `"beginning"`

  :   Set the baseline to 0.

- loess_span:

  Numeric; LOESS span (0, 1\] passed to
  [`predict_baseline_trend`](https://corneliusfritz.github.io/redeem/reference/predict_baseline_trend.md)
  when `baseline_method = "trend"`. Defaults to 0.75.

## Value

A `ranking_redeem` data frame with columns:

- `Cutpoint`:

  Integer value from 0 to `k_max`.

- `Recall`:

  The proportion of test events where the true dyad is ranked at or
  within the cutpoint.

- `Precision`:

  The precision value at the cutpoint.

Additionally, the returned object has the following attributes:

- `"mrr"`:

  Mean Reciprocal Rank (MRR) of the true dyads.

- `"mean_rank"`:

  Mean rank of the true dyads (excluding ranks \> `k_max`).

- `"median_rank"`:

  Median rank of the true dyads (excluding ranks \> `k_max`).

- `"hits_summary"`:

  A data frame summarizing Recall, Precision, and F1 values at K = 1, 5,
  10, and 50.

## Details

For each event observed in the test period (`edgelist_test`), the
function:

1.  Determines the set of all potential candidate dyads (the risk set)
    at that event's timestamp.

2.  Computes the predicted event intensities (or probabilities) for all
    candidate dyads using the fitted model's parameters and the network
    history up to that moment.

3.  Ranks all candidate dyads in descending order of their predicted
    intensities.

4.  Determines the rank of the actually observed dyad.

A well-fitting model will consistently assign higher intensities to the
dyads that actually interact, ranking them near the top.

The function summarizes the rankings across all test events to compute:

- **Mean Reciprocal Rank (MRR)**: The average of the reciprocal ranks of
  the true dyads.

- **Recall at K**: The proportion of test events where the true dyad is
  ranked within the top \\K\\ candidate dyads.

- **Precision at K**: The proportion of top \\K\\ recommendations that
  correspond to true events.
