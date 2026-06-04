# The dem Object

An object of class `dem` returned by the
[`dem`](https://corneliusfritz.github.io/redeem/reference/dem.md)
function, representing a fitted Durational Event Model.

## Value

A `dem` object is a list containing the following components:

- `call`: The matched call.

- `event_numbers`: A vector containing the number of events.

- `model_0_1`: The fitted model for transition 0 -\> 1 (formation).

- `model_1_0`: The fitted model for transition 1 -\> 0 (dissolution).

- `events`: The preprocessed event matrix.

- `formula_0_1`: The formula for transition 0 -\> 1.

- `formula_1_0`: The formula for transition 1 -\> 0.

- `n_nodes`: The number of nodes.

- `simultaneous_interactions`: Logical indicating whether simultaneous
  interactions were allowed.

- `directed`: Logical indicating whether the events are directed.

- `training_start`: The start time of the training period.

- `build_time`: The time at which the estimation dataset started
  building.

- `max_time`: The maximum event time.

- `exogenous_end`: The end time of the exogenous period.

- `time_changepoints`: Time points where baseline intensity changes.

- `labels_changepoints`: Labels for the time intervals.

- `subsample`: Subsample proportion used.

- `return_data`: Logical indicating whether preprocessed data frames
  were returned.

- `runtime`: The estimation runtime.

- `window_map`: The window map used for calculation.

- `preprocessed`: Preprocessed data structures.

## Methods (S3)

The following S3 methods are implemented for `dem` objects:

- `print(x, ...)`: Prints a brief summary of the DEM model.

  - `x`: A `dem` object.

  - `...`: Additional arguments passed to printing function.

- `summary(object, ...)`: Summarizes model results, including parameter
  estimates, standard errors, and fit statistics (AIC/BIC). Returns an
  object of class `summary.dem`.

  - `object`: A `dem` object.

  - `...`: Additional arguments passed to summary method.

- `plot(x, which = 3, separate = FALSE, baseline = FALSE, ...)`:
  Generates trace plots for coefficients and log-likelihood histories.

  - `x`: A `dem` object.

  - `which`: Integer indicating transition plots to display. Options are
    `1` for formation (0 -\> 1), `2` for dissolution (1 -\> 0), or `3`
    (default) for both.

  - `separate`: Logical. If `TRUE`, each plot is generated in a new
    window or as a separate sequence. Defaults to `FALSE`.

  - `baseline`: Logical. If `TRUE`, plots the estimated baseline step
    function. Defaults to `FALSE`.

  - `...`: Additional arguments passed to the underlying trace plotting
    method. Supported parameters include:

    - `coefs`: Logical. If `TRUE` (default), plots iteration traces for
      core/fixed coefficients.

    - `degree`: Logical. If `TRUE` (default), plots iteration traces for
      degree/actor effects.

    - `time`: Logical. If `TRUE` (default), plots iteration traces for
      temporal baseline effects.

    - `llh`: Logical. If `TRUE` (default), plots trace of log-likelihood
      history.

    - `sub_label`: Character. Optional subtitle to display at the bottom
      of the figure.

- `plot_baseline(x, process = c("formation", "dissolution"), ...)`:
  Plots the step function of the estimated baseline intensity.

  - `x`: A `dem` object.

  - `process`: Character. Specifies whether to plot the baseline for
    `"formation"` (0 -\> 1) (default) or the `"dissolution"` (1 -\> 0)
    process.

  - `...`: Additional graphical parameters passed to
    [`plot`](https://rdrr.io/r/graphics/plot.default.html).

- `predict(object, time = NULL, type = c("response", "lp", "terms"), process = c("both", "formation", "dissolution"), ...)`:
  Predicts the intensity, linear predictor, or term contributions for a
  fitted DEM model.

  - `object`: A `dem` object.

  - `time`: Numeric vector; optional time point(s) at which to predict.
    Defaults to NULL.

  - `type`: Character; the type of prediction. Defaults to `"response"`.

  - `process`: Character; the transition process to predict. Defaults to
    `"both"`.

  - `...`: Additional arguments.
