# The rem Object

An object of class `rem` returned by the
[`rem`](https://corneliusfritz.github.io/redeem/reference/rem.md)
function, representing a fitted Relational Event Model.

## Value

A `rem` object is a list containing the following components:

- `call`: The matched call.

- `event_numbers`: A vector containing the number of events.

- `model`: The fitted underlying model.

- `events`: The preprocessed event matrix.

- `formula`: The formula used.

- `n_nodes`: The number of nodes.

- `directed`: Logical indicating whether the events are directed.

- `build_time`: The time at which the estimation dataset started
  building.

- `max_time`: The maximum event time.

- `time_changepoints`: Time points where baseline intensity changes.

- `labels_changepoints`: Labels for the time intervals.

- `training_start`: The start time of the training period.

- `exogenous_end`: The end time of the exogenous period.

- `subsample`: Subsample proportion used.

- `return_data`: Logical indicating whether preprocessed data frames
  were returned.

- `runtime`: The estimation runtime.

- `window_map`: The window map used for calculation.

- `preprocessed`: Preprocessed data structures.

## Methods (S3)

The following S3 methods are implemented for `rem` objects:

- `print(x, ...)`: Prints a brief summary of the fitted REM model.

  - `x`: A `rem` object.

  - `...`: Additional arguments passed to printing function.

- `summary(object, ...)`: Summarizes model results, including parameter
  estimates, standard errors, and fit statistics. Returns an object of
  class `summary.rem`.

  - `object`: A `rem` object.

  - `...`: Additional arguments passed to summary method.

- `plot(x, baseline = FALSE, ...)`: Generates trace plots for model
  coefficients and log-likelihood histories.

  - `x`: A `rem` object.

  - `baseline`: Logical. If `TRUE`, plots the estimated baseline step
    function. Defaults to `FALSE`.

  - `...`: Additional arguments passed to underlying trace plotting
    method. Supported parameters include:

    - `coefs`: Logical. If `TRUE` (default), plots iteration traces for
      core/fixed coefficients.

    - `degree`: Logical. If `TRUE` (default), plots iteration traces for
      degree/actor effects.

    - `time`: Logical. If `TRUE` (default), plots iteration traces for
      temporal baseline effects.

    - `llh`: Logical. If `TRUE` (default), plots trace of log-likelihood
      history.

    - `separate`: Logical. If `TRUE`, each plot is generated in a new
      window or as a separate sequence. Defaults to `FALSE`.

    - `sub_label`: Character. Optional subtitle to display at the bottom
      of the figure.

- `logLik(object, ...)`: Extracts the log-likelihood of the fitted
  model.

  - `object`: A `rem` object.

  - `...`: Additional arguments.

- `plot_baseline(x, ...)`: Plots the step function of the estimated
  baseline intensity.

  - `x`: A `rem` object.

  - `...`: Additional graphical parameters passed to
    [`plot`](https://rdrr.io/r/graphics/plot.default.html).

- `predict(object, time = NULL, type = c("response", "lp", "terms"), ...)`:
  Predicts the intensity, linear predictor, or term contributions for a
  fitted REM model.

  - `object`: A `rem` object.

  - `time`: Numeric vector; optional time point(s) at which to predict.
    Defaults to NULL.

  - `type`: Character; the type of prediction. Defaults to `"response"`.

  - `...`: Additional arguments.
