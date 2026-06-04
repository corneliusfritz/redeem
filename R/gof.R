utils::globalVariables(c("prediction", "event", "pair_id", ".SD", ".", "from", "to", "status", "time_cat", "effect", "time_new", "grid_start", "grid_end"))

#' Predict the baseline intensity trend at one or more future time points
#'
#' Decomposes the estimated piecewise-constant log-baseline (\code{est_time})
#' into a smooth trend component (via LOESS) and a seasonal/residual component,
#' following the same approach used in the application plot script.  The fitted
#' trend is then \emph{extrapolated} to \code{target_times} using
#' \code{predict.loess()} so that the baseline used for out-of-sample scoring
#' reflects the long-run level of activity rather than any arbitrary fixed
#' value.
#'
#' @param model A \code{redeem_result} object with non-null \code{est_time}
#'   and \code{time_changepoints} fields.
#' @param target_times Numeric vector: the times at which to predict the trend
#'   (typically the unique timestamps of the test events).
#' @param loess_span Numeric; the span argument passed to \code{stats::loess}.
#'   Larger values give a smoother (more conservative) trend extrapolation. Defaults to 0.75.
#'
#' @return A numeric vector of predicted log-baseline (trend component) values,
#'   one per element of \code{target_times}.  Falls back to
#'   \code{mean(est_time)} for each time point if there are fewer than 3
#'   observations or if the LOESS fit fails.
#'
#' @details
#' The decomposition mirrors the plot code in the application:
#' \enumerate{
#'   \item Build a data frame of \code{(time, est_time)} using the changepoints
#'         stored in \code{model$time_changepoints}.  The first interval [0,
#'         changepoint_1) is given time = 0; each subsequent interval gets the
#'         corresponding changepoint value.
#'   \item Fit LOESS on the log-scale \code{est_time} values.
#'   \item Predict at each \code{target_times}; predictions are clamped to the
#'         range of the observed \code{est_time} to avoid wild extrapolation.
#' }
#' @export
predict_baseline_trend <- function(model, target_times, loess_span = 0.75) {
  # Guard: need est_time
  est_time <- model$est_time
  if (is.null(est_time) || length(est_time) == 0) {
    return(rep(0.0, length(target_times)))
  }

  full_baseline <- isTRUE(model$full_baseline)

  # Build (time, value) pairs on the log scale (est_time is already log-scale)
  # For full_baseline: est_time has K+1 entries aligned with [0, cp1, cp2, ...]
  # For normal:        est_time has K entries; a 0 is prepended for the [0, cp1) bin
  gamma <- if (full_baseline) est_time else c(0, est_time)

  changepoints <- model$time_changepoints
  # Corresponding representative times for each interval
  # Interval 1 -> time 0, Interval k+1 -> changepoints[k]
  times <- c(0, changepoints)

  # Align lengths (safety)
  n <- min(length(gamma), length(times))
  fallback <- mean(est_time, na.rm = TRUE)
  if (n < 3) {
    # Not enough points for trend - fall back to mean for all target times
    return(rep(fallback, length(target_times)))
  }
  times <- times[seq_len(n)]
  gamma <- gamma[seq_len(n)]

  # Filter out non-finite values
  ok <- is.finite(gamma) & is.finite(times)
  if (sum(ok) < 3) {
    return(rep(mean(est_time[is.finite(est_time)], na.rm = TRUE), length(target_times)))
  }
  times_ok <- times[ok]
  gamma_ok <- gamma[ok]

  # 1. Fit LOESS trend
  trend_fit <- tryCatch(
    stats::loess(gamma_ok ~ times_ok,
      span = loess_span,
      control = stats::loess.control(surface = "direct")
    ),
    error = function(e) NULL
  )
  if (is.null(trend_fit)) {
    return(rep(mean(gamma_ok), length(target_times)))
  }

  # Predict trend on the training points
  trend_fitted <- tryCatch(
    stats::predict(trend_fit, newdata = data.frame(times_ok = times_ok)),
    error = function(e) rep(mean(gamma_ok), length(times_ok))
  )
  trend_fitted[is.na(trend_fitted) | !is.finite(trend_fitted)] <- mean(gamma_ok)

  # Compute residuals for seasonal estimation
  res_ok <- gamma_ok - trend_fitted

  # 2. Fit Daily and Weekly seasonal patterns on residuals using Harmonic Regression
  # Period constants in seconds
  period_day <- 86400
  period_week <- 604800

  time_range <- diff(range(times_ok))
  fit_seasonal <- FALSE
  use_weekly <- FALSE

  # Build the seasonal model formulas dynamically based on dataset span
  if (time_range >= period_week) {
    fit_seasonal <- TRUE
    use_weekly <- TRUE
  } else if (time_range >= period_day) {
    fit_seasonal <- TRUE
    use_weekly <- FALSE
  }

  seasonal_fit <- NULL
  if (fit_seasonal) {
    t_day <- times_ok %% period_day
    df_fit <- data.frame(
      res = res_ok,
      sin_d1 = sin(2 * pi * t_day / period_day),
      cos_d1 = cos(2 * pi * t_day / period_day),
      sin_d2 = sin(4 * pi * t_day / period_day),
      cos_d2 = cos(4 * pi * t_day / period_day)
    )
    formula_str <- "res ~ sin_d1 + cos_d1 + sin_d2 + cos_d2"

    if (use_weekly) {
      t_week <- times_ok %% period_week
      df_fit$sin_w1 <- sin(2 * pi * t_week / period_week)
      df_fit$cos_w1 <- cos(2 * pi * t_week / period_week)
      df_fit$sin_w2 <- sin(4 * pi * t_week / period_week)
      df_fit$cos_w2 <- cos(4 * pi * t_week / period_week)
      formula_str <- paste(formula_str, "+ sin_w1 + cos_w1 + sin_w2 + cos_w2")
    }

    seasonal_fit <- tryCatch(
      stats::lm(stats::as.formula(formula_str), data = df_fit),
      error = function(e) NULL
    )
  }

  # 3. Predict Trend at target_times
  pred_trend <- tryCatch(
    stats::predict(trend_fit, newdata = data.frame(times_ok = target_times)),
    error = function(e) rep(NA_real_, length(target_times)),
    warning = function(w) {
      suppressWarnings(
        stats::predict(trend_fit, newdata = data.frame(times_ok = target_times))
      )
    }
  )

  # Clamp each predicted trend to the observed trend range to avoid wild extrapolation
  obs_range <- range(gamma_ok)
  pred_trend <- ifelse(
    is.na(pred_trend) | !is.finite(pred_trend),
    mean(gamma_ok),
    pmax(obs_range[1], pmin(obs_range[2], pred_trend))
  )

  # 4. Predict Seasonality at target_times
  pred_seasonal <- rep(0.0, length(target_times))
  if (!is.null(seasonal_fit)) {
    t_day_target <- target_times %% period_day
    df_pred <- data.frame(
      sin_d1 = sin(2 * pi * t_day_target / period_day),
      cos_d1 = cos(2 * pi * t_day_target / period_day),
      sin_d2 = sin(4 * pi * t_day_target / period_day),
      cos_d2 = cos(4 * pi * t_day_target / period_day)
    )
    if (use_weekly) {
      t_week_target <- target_times %% period_week
      df_pred$sin_w1 <- sin(2 * pi * t_week_target / period_week)
      df_pred$cos_w1 <- cos(2 * pi * t_week_target / period_week)
      df_pred$sin_w2 <- sin(4 * pi * t_week_target / period_week)
      df_pred$cos_w2 <- cos(4 * pi * t_week_target / period_week)
    }

    pred_seasonal <- tryCatch(
      stats::predict(seasonal_fit, newdata = df_pred),
      error = function(e) rep(0.0, length(target_times))
    )
    pred_seasonal[is.na(pred_seasonal) | !is.finite(pred_seasonal)] <- 0.0
  }

  # 5. Combine Trend and Seasonality
  predicted <- pred_trend + pred_seasonal

  return(as.double(predicted))
}

.make_baseline_vec <- function(model_obj, unique_times, method, loess_span = 0.75) {
  n <- length(unique_times)
  if (is.null(model_obj) || is.null(model_obj$est_time) ||
    length(model_obj$est_time) == 0) {
    return(rep(0.0, n))
  }
  switch(method,
    "trend" = predict_baseline_trend(model_obj,
      target_times = unique_times,
      loess_span = loess_span
    ),
    "mean" = rep(mean(model_obj$est_time, na.rm = TRUE), n),
    "beginning" = rep(0, n),
    "last" = {
      gamma <- if (isTRUE(model_obj$full_baseline)) {
        model_obj$est_time
      } else {
        c(0, model_obj$est_time)
      }
      rep(tail(gamma[is.finite(gamma)], 1), n)
    }
  )
}


#' Get ranking for test events (Out-of-Sample Goodness-of-Fit)
#'
#' @description
#' Evaluates the out-of-sample predictive performance of a fitted model on a test event sequence
#' using a ranking-based Goodness-of-Fit (GoF) procedure.
#'
#' @details
#' For each event observed in the test period (\code{edgelist_test}), the function:
#' \enumerate{
#'   \item Determines the set of all potential candidate dyads (the risk set) at that event's timestamp.
#'   \item Computes the predicted event intensities (or probabilities) for all candidate dyads using the fitted model's parameters and the network history up to that moment.
#'   \item Ranks all candidate dyads in descending order of their predicted intensities.
#'   \item Determines the rank of the actually observed dyad.
#' }
#' A well-fitting model will consistently assign higher intensities to the dyads that actually interact, ranking them near the top.
#'
#' The function summarizes the rankings across all test events to compute:
#' \itemize{
#'   \item \strong{Mean Reciprocal Rank (MRR)}: The average of the reciprocal ranks of the true dyads.
#'   \item \strong{Recall at K}: The proportion of test events where the true dyad is ranked within the top \eqn{K} candidate dyads.
#'   \item \strong{Precision at K}: The proportion of top \eqn{K} recommendations that correspond to true events.
#' }
#'
#' @param object A \code{redeem} object (either \code{\link{rem}} or \code{\link{dem}}).
#' @param verbose Logical; if `TRUE`, prints verbose output. Defaults to FALSE.
#' @param k_max Maximum number of ranked pairs to return. Defaults to 1000.
#' @param edgelist_test A matrix of test events (timing, from, to, type).
#' @param edgelist_train A matrix of train events (timing, from, to, type). Defaults to NULL.
#' @param baseline_method Character; how to compute the fixed log-baseline
#'   intensity used for out-of-sample scoring. Defaults to \code{"trend"}. One of:
#'   \describe{
#'     \item{\code{"trend"}}{Fit a LOESS trend to the estimated
#'       piecewise-constant log-baseline (\code{est_time}) over training time
#'       and extrapolate it to the start of the test period.  This mirrors the
#'       trend decomposition used in the application plot script and typically
#'       yields a better forecast than a fixed mean.}
#'     \item{\code{"mean"}}{Use the simple mean of \code{est_time}.}
#'     \item{\code{"last"}}{Use the last estimated log-baseline value
#'       (i.e.\ the value from the most recent training interval).}
#'     \item{\code{"beginning"}}{Set the baseline to 0.}
#'   }
#' @param loess_span Numeric; LOESS span (0, 1] passed to
#'   \code{\link{predict_baseline_trend}} when \code{baseline_method = "trend"}. Defaults to 0.75.
#' @param ties.method Character; the method to handle ties when ranking event intensities,
#'   passed directly to \code{\link[base]{rank}}. Defaults to \code{"average"}. One of:
#'   \describe{
#'     \item{\code{"average"}}{Assigns the average of the ranks of all tied elements to each.}
#'     \item{\code{"first"}}{Breaks ties by the order they appear in the data structure.}
#'     \item{\code{"last"}}{Breaks ties by the reverse order of their appearance.}
#'     \item{\code{"random"}}{Breaks ties randomly, ensuring no systematic bias.}
#'     \item{\code{"max"}}{Assigns the maximum of the ranks of the tied elements to all.}
#'     \item{\code{"min"}}{Assigns the minimum of the ranks of the tied elements to all.}
#'   }
#' @param return_probabilities Logical; if TRUE, returns the predicted probabilities/scores instead of recall curves. Defaults to FALSE.
#'
#' @return A \code{ranking_redeem} data frame with columns:
#'   \describe{
#'     \item{\code{Cutpoint}}{Integer value from 0 to \code{k_max}.}
#'     \item{\code{Recall}}{The proportion of test events where the true dyad is ranked at or within the cutpoint.}
#'     \item{\code{Precision}}{The precision value at the cutpoint.}
#'   }
#'   Additionally, the returned object has the following attributes:
#'   \describe{
#'     \item{\code{"mrr"}}{Mean Reciprocal Rank (MRR) of the true dyads.}
#'     \item{\code{"mean_rank"}}{Mean rank of the true dyads (excluding ranks > \code{k_max}).}
#'     \item{\code{"median_rank"}}{Median rank of the true dyads (excluding ranks > \code{k_max}).}
#'     \item{\code{"hits_summary"}}{A data frame summarizing Recall, Precision, and F1 values at K = 1, 5, 10, and 50.}
#'   }
#' @export
get_ranking <- function(object,
                        verbose = FALSE,
                        k_max = 1000,
                        edgelist_test,
                        edgelist_train = NULL,
                        ties.method = c("average", "first", "last", "random", "max", "min"),
                        return_probabilities = FALSE,
                        baseline_method = c("trend", "mean", "last", "beginning"),
                        loess_span = 0.75) {
  ties.method <- match.arg(ties.method)
  baseline_method <- match.arg(baseline_method)

  if (is.null(edgelist_test) || nrow(edgelist_test) == 0) {
    res <- data.frame(
      Cutpoint = 0:k_max,
      Recall = rep(0, k_max + 1),
      Precision = rep(0, k_max + 1)
    )
    class(res) <- c("ranking_redeem", class(res))
    attr(res, "mrr") <- 0
    attr(res, "mean_rank") <- NA
    attr(res, "median_rank") <- NA
    attr(res, "hits_summary") <- data.frame(
      Metric = character(0),
      Value = numeric(0),
      stringsAsFactors = FALSE
    )
    return(res)
  }


  formula_0_1 <- if (inherits(object, "rem")) object$formula else object$formula_0_1
  formula_1_0 <- if (inherits(object, "rem")) NULL else object$formula_1_0
  n_nodes <- object$n_nodes

  model_0_1 <- if (inherits(object, "rem")) object$model else object$model_0_1
  model_1_0 <- if (inherits(object, "rem")) NULL else object$model_1_0

  if (is.null(model_0_1)) {
    stop("model_0_1 is NULL. Check if object is a valid rem or dem result.")
  }
  coef_0_1 <- if (!is.null(model_0_1$est_core)) model_0_1$est_core else numeric(0)
  coef_1_0 <- if (!is.null(model_1_0) && !is.null(model_1_0$est_core)) model_1_0$est_core else numeric(0)

  coef_0_1_degree <- if (!is.null(model_0_1$est_degree)) model_0_1$est_degree else numeric(0)
  coef_1_0_degree <- if (!is.null(model_1_0) && !is.null(model_1_0$est_degree)) model_1_0$est_degree else numeric(0)

  # Prevent out-of-sample -Inf penalties by setting -Inf values to the minimum of finite values
  if (length(coef_0_1_degree) > 0) {
    finite_deg <- coef_0_1_degree[is.finite(coef_0_1_degree)]
    min_finite <- if (length(finite_deg) > 0) min(finite_deg) else 0.0
    coef_0_1_degree[is.infinite(coef_0_1_degree) & coef_0_1_degree < 0] <- min_finite
  }
  if (length(coef_1_0_degree) > 0) {
    finite_deg <- coef_1_0_degree[is.finite(coef_1_0_degree)]
    min_finite <- if (length(finite_deg) > 0) min(finite_deg) else 0.0
    coef_1_0_degree[is.infinite(coef_1_0_degree) & coef_1_0_degree < 0] <- min_finite
  }


  # Guard against semiparametric models
  if (inherits(model_0_1, "coxph") || inherits(model_1_0, "coxph") ||
    isTRUE(object$semiparametric) || inherits(object, "dem.cox")) {
    stop("get_ranking is currently not implemented for semiparametric models.")
  }

  simultaneous_interactions <- if (is.null(object$simultaneous_interactions)) FALSE else object$simultaneous_interactions

  preprocessed <- if (!is.null(object[["preprocessed"]])) {
    object[["preprocessed"]]
  } else {
    formula_preprocess(
      formula_1_0 = formula_1_0,
      formula_0_1 = formula_0_1,
      events = object$events,
      n_nodes = n_nodes,
      model_type = if (inherits(object, "rem")) "rem" else "dem",
      directed = object$directed
    )
  }

  if (is.null(edgelist_train)) {
    if (!is.null(object$events)) {
      edgelist_train <- object$events
    } else if (!is.null(preprocessed$events)) {
      edgelist_train <- preprocessed$events
    } else {
      stop("No edgelist_train or events found in the object.")
    }
  }

  edgelist_test <- as.matrix(edgelist_test)
  edgelist_train <- as.matrix(edgelist_train)

  if (ncol(edgelist_test) == 3) {
    edgelist_test <- cbind(edgelist_test, 1)
  }
  if (ncol(edgelist_train) == 3) {
    edgelist_train <- cbind(edgelist_train, 1)
  }

  if (!object$directed) {
    if (nrow(edgelist_test) > 0) {
      swap_test <- edgelist_test[, 2] > edgelist_test[, 3]
      if (any(swap_test)) {
        tmp_vals <- edgelist_test[swap_test, 2]
        edgelist_test[swap_test, 2] <- edgelist_test[swap_test, 3]
        edgelist_test[swap_test, 3] <- tmp_vals
      }
    }
    if (nrow(edgelist_train) > 0) {
      swap_train <- edgelist_train[, 2] > edgelist_train[, 3]
      if (any(swap_train)) {
        tmp_vals <- edgelist_train[swap_train, 2]
        edgelist_train[swap_train, 2] <- edgelist_train[swap_train, 3]
        edgelist_train[swap_train, 3] <- tmp_vals
      }
    }
  }

  # Check for multi-stream
  is_multistream <- any(!sapply(preprocessed$stream_list, is.null))
  if (is_multistream) {
    stop("get_ranking() does not yet support multi-stream event models.")
  }

  coef_1_0_augmented <- numeric(length = length(preprocessed$coef_names))
  idx_1_0 <- match(preprocessed$preprocess_1_0$coef_names, preprocessed$coef_names)
  if (length(idx_1_0) > 0 && length(coef_1_0) > 0) {
    coef_1_0_augmented[idx_1_0] <- coef_1_0
  }

  coef_0_1_augmented <- numeric(length = length(preprocessed$coef_names))
  idx_0_1 <- match(preprocessed$preprocess_0_1$coef_names, preprocessed$coef_names)
  if (length(idx_0_1) > 0 && length(coef_0_1) > 0) {
    coef_0_1_augmented[idx_0_1] <- coef_0_1
  }

  unique_test_times <- sort(unique(edgelist_test[, 1]))

  # plot(predict_baseline_trend(model_0_1,
  #                        target_times = unique_test_times,
  #                        loess_span = 0.1))

  baseline_0_1_vec <- .make_baseline_vec(model_0_1, unique_test_times, baseline_method, loess_span)
  baseline_1_0_vec <- .make_baseline_vec(model_1_0, unique_test_times, baseline_method, loess_span)

  # Call C++ helper for ranking
  tmp <- get_probabilities_per_test_event(
    terms = as.character(unlist(preprocessed$term_names)),
    data_list = lapply(preprocessed$data_list, function(x) {
      if (is.list(x)) {
        return(x)
      } # Preserve list for TV
      m <- as.matrix(x)
      storage.mode(m) <- "double"
      m
    }),
    transformations = as.character(unlist(preprocessed$transformation_list)),
    n_nodes = as.integer(preprocessed$n_nodes),
    verbose = as.logical(verbose),
    directed = as.logical(object$directed),
    coef_0_1 = as.double(coef_0_1_augmented),
    coef_1_0 = as.double(coef_1_0_augmented),
    degree_coef_0_1 = as.double(coef_0_1_degree),
    degree_coef_1_0 = as.double(coef_1_0_degree),
    simultaneous_interactions = as.logical(simultaneous_interactions),
    edgelist_train = as.matrix(edgelist_train),
    edgelist_test = as.matrix(edgelist_test),
    k = as.integer(k_max),
    is_rem = as.logical(inherits(object, "rem")),
    window_info = if (length(preprocessed$window_map) > 0) {
      as.list(stats::setNames(as.numeric(names(preprocessed$window_map)), preprocessed$window_map))
    } else {
      list()
    },
    baseline_0_1 = as.double(baseline_0_1_vec),
    baseline_1_0 = as.double(baseline_1_0_vec)
  )


  if (isTRUE(return_probabilities)) {
    return(tmp)
  }

  info <- unlist(lapply(tmp, function(x) {
    keys_observed <- paste(x$observed[, 1], x$observed[, 2], sep = "_")
    keys_predicted <- paste(x$predicted[, 1], x$predicted[, 2], sep = "_")
    idx <- match(keys_observed, keys_predicted)
    if (ties.method == "first") {
      return(idx)
    } else {
      ranks <- rank(-x$predicted[, 3], ties.method = ties.method)
      return(ranks[idx])
    }
  }))
  info[is.na(info)] <- k_max + 1

  # Compute additional ranking metrics from resolved ranks
  rr <- ifelse(info <= k_max, 1 / info, 0)
  mrr <- mean(rr, na.rm = TRUE)

  recalled_ranks <- info[info <= k_max]
  mean_rank <- if (length(recalled_ranks) > 0) mean(recalled_ranks) else NA
  median_rank <- if (length(recalled_ranks) > 0) stats::median(recalled_ranks) else NA

  hits_at <- c(1, 5, 10, 50)
  hits_at <- hits_at[hits_at <= k_max]

  recall_vals <- sapply(hits_at, function(k) {
    mean(info <= k)
  })
  names(recall_vals) <- paste0("Recall@", hits_at)

  precision_vals <- sapply(hits_at, function(k) {
    mean(info <= k) / k
  })
  names(precision_vals) <- paste0("Precision@", hits_at)

  f1_vals <- sapply(hits_at, function(k) {
    rec <- mean(info <= k)
    prec <- rec / k
    if (rec + prec > 0) {
      2 * (prec * rec) / (prec + rec)
    } else {
      0
    }
  })
  names(f1_vals) <- paste0("F1@", hits_at)

  hits_summary <- data.frame(
    Metric = c(paste0("Recall@", hits_at), paste0("Precision@", hits_at), paste0("F1@", hits_at)),
    Value = c(recall_vals, precision_vals, f1_vals),
    stringsAsFactors = FALSE
  )

  res <- data.frame(
    Cutpoint = 0:k_max,
    Recall = c(0, findInterval(x = 1:k_max, sort(info)) / length(sort(info)))
  )
  res$Precision <- c(0, res$Recall[-1] / (1:k_max))
  class(res) <- c("ranking_redeem", class(res))
  attr(res, "mrr") <- mrr
  attr(res, "mean_rank") <- mean_rank
  attr(res, "median_rank") <- median_rank
  attr(res, "hits_summary") <- hits_summary
  return(res)
}


#' Out-of-sample Log-Likelihood (Proper Scoring Rule)
#'
#' This function computes the out-of-sample log-likelihood (a strictly proper scoring rule) for each test event
#' under a fitted REM or DEM.
#'
#' @param object A \code{redeem} object (either \code{\link{rem}} or \code{\link{dem}}).
#' @param verbose Logical; if `TRUE`, prints verbose output. Defaults to FALSE.
#' @param edgelist_test A matrix or data frame of test events (timing, from, to, type).
#' @param edgelist_train A matrix or data frame of train events (timing, from, to, type). Defaults to `NULL`,
#'   in which case it retrieves the training events from the `object` or the preprocessed data.
#' @param baseline_method Character; how to compute the fixed log-baseline
#'   intensity used for out-of-sample scoring. One of: `"last"` (uses
#'   the last estimated baseline value), `"trend"` (extrapolates a LOESS trend),
#'   `"mean"`, or `"beginning"`. Defaults to `"last"`.
#' @param loess_span Numeric; LOESS span (0, 1] passed to
#'   \code{\link{predict_baseline_trend}} when \code{baseline_method = "trend"}. Defaults to 0.75.
#'
#' @seealso \code{\link{rem_object}} and \code{\link{dem_object}} for details on prediction methods.
#'
#' @return A numeric vector of log-likelihoods for each test event.
#' @export
get_oos_likelihood <- function(object,
                               verbose = FALSE,
                               edgelist_test,
                               edgelist_train = NULL,
                               baseline_method = c("last", "trend", "mean", "beginning"),
                               loess_span = 0.75) {
  baseline_method <- match.arg(baseline_method)
  formula_0_1 <- if (inherits(object, "rem")) object$formula else object$formula_0_1
  formula_1_0 <- if (inherits(object, "rem")) NULL else object$formula_1_0
  n_nodes <- object$n_nodes

  preprocessed <- if (!is.null(object[["preprocessed"]])) {
    object[["preprocessed"]]
  } else {
    formula_preprocess(
      formula_1_0 = formula_1_0,
      formula_0_1 = formula_0_1,
      events = object$events,
      n_nodes = n_nodes,
      model_type = if (inherits(object, "rem")) "rem" else "dem",
      directed = object$directed
    )
  }

  if (is.null(edgelist_train)) {
    if (!is.null(object$events)) {
      edgelist_train <- object$events
    } else if (!is.null(preprocessed$events)) {
      edgelist_train <- preprocessed$events
    } else {
      stop("No edgelist_train or events found in the object.")
    }
  }

  edgelist_test <- as.matrix(edgelist_test)
  edgelist_train <- as.matrix(edgelist_train)

  if (ncol(edgelist_test) == 3) {
    edgelist_test <- cbind(edgelist_test, 1)
  }
  if (ncol(edgelist_train) == 3) {
    edgelist_train <- cbind(edgelist_train, 1)
  }

  if (!object$directed) {
    if (nrow(edgelist_test) > 0) {
      swap_test <- edgelist_test[, 2] > edgelist_test[, 3]
      if (any(swap_test)) {
        tmp_vals <- edgelist_test[swap_test, 2]
        edgelist_test[swap_test, 2] <- edgelist_test[swap_test, 3]
        edgelist_test[swap_test, 3] <- tmp_vals
      }
    }
    if (nrow(edgelist_train) > 0) {
      swap_train <- edgelist_train[, 2] > edgelist_train[, 3]
      if (any(swap_train)) {
        tmp_vals <- edgelist_train[swap_train, 2]
        edgelist_train[swap_train, 2] <- edgelist_train[swap_train, 3]
        edgelist_train[swap_train, 3] <- tmp_vals
      }
    }
  }

  model_0_1 <- if (inherits(object, "rem")) object$model else object$model_0_1
  model_1_0 <- if (inherits(object, "rem")) NULL else object$model_1_0

  if (is.null(model_0_1)) {
    stop("model_0_1 is NULL. Check if object is a valid rem or dem result.")
  }
  coef_0_1 <- if (!is.null(model_0_1$est_core)) model_0_1$est_core else numeric(0)
  coef_1_0 <- if (!is.null(model_1_0) && !is.null(model_1_0$est_core)) model_1_0$est_core else numeric(0)

  coef_0_1_degree <- if (!is.null(model_0_1$est_degree)) model_0_1$est_degree else numeric(0)
  coef_1_0_degree <- if (!is.null(model_1_0) && !is.null(model_1_0$est_degree)) model_1_0$est_degree else numeric(0)

  # Prevent out-of-sample -Inf penalties by setting -Inf values to the minimum of finite values
  if (length(coef_0_1_degree) > 0) {
    finite_deg <- coef_0_1_degree[is.finite(coef_0_1_degree)]
    min_finite <- if (length(finite_deg) > 0) min(finite_deg) else 0.0
    coef_0_1_degree[is.infinite(coef_0_1_degree) & coef_0_1_degree < 0] <- min_finite
  }
  if (length(coef_1_0_degree) > 0) {
    finite_deg <- coef_1_0_degree[is.finite(coef_1_0_degree)]
    min_finite <- if (length(finite_deg) > 0) min(finite_deg) else 0.0
    coef_1_0_degree[is.infinite(coef_1_0_degree) & coef_1_0_degree < 0] <- min_finite
  }

  # Guard against semiparametric models
  if (inherits(model_0_1, "coxph") || inherits(model_1_0, "coxph") ||
    isTRUE(object$semiparametric) || inherits(object, "dem.cox")) {
    stop("get_oos_likelihood is currently not implemented for semiparametric models.")
  }

  simultaneous_interactions <- if (is.null(object$simultaneous_interactions)) FALSE else object$simultaneous_interactions

  # Check for multi-stream
  is_multistream <- any(!sapply(preprocessed$stream_list, is.null))
  if (is_multistream) {
    stop("get_oos_likelihood() does not yet support multi-stream event models.")
  }

  coef_0_1_augmented <- numeric(length = length(preprocessed$coef_names))
  idx_0_1 <- match(preprocessed$preprocess_0_1$coef_names, preprocessed$coef_names)
  if (length(idx_0_1) > 0 && length(coef_0_1) > 0) {
    coef_0_1_augmented[idx_0_1] <- coef_0_1
  }

  coef_1_0_augmented <- numeric(length = length(preprocessed$coef_names))
  idx_1_0 <- match(preprocessed$preprocess_1_0$coef_names, preprocessed$coef_names)
  if (length(idx_1_0) > 0 && length(coef_1_0) > 0) {
    coef_1_0_augmented[idx_1_0] <- coef_1_0
  }

  # Compute per-unique-test-time baseline vectors (same logic as get_ranking)
  unique_test_times <- sort(unique(edgelist_test[, 1]))

  baseline_0_1_vec <- .make_baseline_vec(model_0_1, unique_test_times, baseline_method, loess_span)
  baseline_1_0_vec <- .make_baseline_vec(model_1_0, unique_test_times, baseline_method, loess_span)

  ll_vals <- get_oos_likelihood_cpp(
    terms = as.character(unlist(preprocessed$term_names)),
    data_list = lapply(preprocessed$data_list, function(x) {
      if (is.list(x)) {
        return(x)
      } # Preserve list for TV
      m <- as.matrix(x)
      storage.mode(m) <- "double"
      m
    }),
    transformations = as.character(unlist(preprocessed$transformation_list)),
    n_nodes = as.integer(preprocessed$n_nodes),
    verbose = as.logical(verbose),
    directed = as.logical(object$directed),
    coef_0_1 = as.double(coef_0_1_augmented),
    coef_1_0 = as.double(coef_1_0_augmented),
    degree_coef_0_1 = as.double(coef_0_1_degree),
    degree_coef_1_0 = as.double(coef_1_0_degree),
    simultaneous_interactions = as.logical(simultaneous_interactions),
    edgelist_train = as.matrix(edgelist_train),
    edgelist_test = as.matrix(edgelist_test),
    is_rem = as.logical(inherits(object, "rem")),
    window_info = if (length(preprocessed$window_map) > 0) {
      as.list(stats::setNames(as.numeric(names(preprocessed$window_map)), preprocessed$window_map))
    } else {
      list()
    },
    baseline_0_1 = as.double(baseline_0_1_vec),
    baseline_1_0 = as.double(baseline_1_0_vec)
  )


  return(ll_vals)
}


#' Get residuals for model diagnostics (Cox-Snell Residuals)
#'
#' @description
#' Computes Cox-Snell residuals for a fitted model to diagnose
#' goodness-of-fit and calibration.
#'
#' @details
#' Cox-Snell residuals are a standard diagnostic tool for continuous-time
#' survival models and counting processes.
#' Under the true model specification, the integrated cumulative intensity
#' computed up to the exact time of an observed event
#' is distributed as a standard exponential random variable, i.e.,
#' \eqn{\Lambda_{ij}(t_k) \sim Exp(1)}.
#'
#' Consequently, if the model is correctly specified:
#' \itemize{
#'   \item The empirical survival function of these residuals should
#'     closely match the theoretical survival function of a standard
#'     exponential distribution, \eqn{S(r) = \exp(-r)}.
#'   \item Deviations between the empirical Kaplan-Meier curve of the
#'     residuals and the theoretical exponential curve signal model
#'     misspecification, unmodeled dyadic heterogeneity, or
#'     non-stationarity.
#' }
#'
#' The function can compute residuals for both the formation/incidence
#' (\eqn{0 \rightarrow 1}) process and the dissolution/duration
#' (\eqn{1 \rightarrow 0}) process.
#'
#' @references
#' Cox, D. R., & Snell, E. J. (1968). A general definition of residuals. Journal of the Royal Statistical Society: Series B (Methodological), 30(2), 248-265.
#'
#' @param object A \code{redeem} object (either \code{\link{rem}} or
#'   \code{\link{dem}}).
#' @param get_0_1 Logical; if `TRUE`, computes residuals for the
#'   formation (0 -> 1) process. Defaults to TRUE.
#' @param get_1_0 Logical; if `TRUE`, computes residuals for the
#'   dissolution (1 -> 0) process. Defaults to TRUE.
#' @param raw Logical; if `TRUE`, returns the raw Cox-Snell residuals. Defaults to FALSE.
#'
#' @return If `raw = TRUE`, a list containing the raw residuals for the
#'   selected process(es).
#'   If `raw = FALSE`, a list of data frames containing the Kaplan-Meier
#'   coordinates (`time`, `surv`) and the corresponding `theoretical`
#'   standard exponential survival values.
#'
#' @export
get_residuals <- function(object, get_0_1 = TRUE, get_1_0 = TRUE, raw = FALSE) {
  model_0_1 <- if (inherits(object, "rem")) object$model else object$model_0_1
  model_1_0 <- if (inherits(object, "rem")) NULL else object$model_1_0


  # Reconstruct data if missing
  reproduced_data <- NULL
  if ((get_0_1 && !is.null(model_0_1) && is.null(model_0_1$data)) ||
    (get_1_0 && !is.null(model_1_0) && is.null(model_1_0$data))) {
    # Guard against semiparametric models
    if ((!is.null(model_0_1) && inherits(model_0_1, "coxph")) ||
      (!is.null(model_1_0) && inherits(model_1_0, "coxph")) ||
      isTRUE(object$semiparametric)) {
      stop("get_residuals is currently not implemented for semiparametric models.")
    }
    reproduced_data <- reproduce_model_data(object)
  }

  resid_data_0_1 <- NULL
  if (get_0_1 && !is.null(model_0_1)) {
    data_0_1 <- if (!is.null(model_0_1$data)) {
      model_0_1$data
    } else if (!is.null(reproduced_data)) {
      if (is.data.frame(reproduced_data)) reproduced_data else reproduced_data$data_0_1
    } else {
      NULL
    }

    if (!is.null(data_0_1)) {
      dt_0_1 <- data.table::as.data.table(data_0_1)

      # Determine column names inside j - grouping by pair_id
      resid_data_0_1 <- dt_0_1[,
        {
          # Access the prediction column if exists, otherwise compute it
          curr_pred_vals <- if ("prediction" %in% names(.SD)) .SD$prediction else calculate_predictions_helper(model_0_1, .SD)

          cum_lambda <- cumsum(as.double(curr_pred_vals))
          event_indices <- which(event == 1)
          final_val <- if (length(cum_lambda) > 0) tail(cum_lambda, 1) else 0

          if (length(event_indices) > 0) {
            event_vals <- cum_lambda[event_indices]
            event_resids <- diff(c(0, event_vals))
            status_vec <- rep(1, length(event_resids))

            tail_resid <- final_val - tail(event_vals, 1)
            if (abs(tail_resid) < 1e-10) {
              list(residual = event_resids, status = status_vec)
            } else {
              list(residual = c(event_resids, tail_resid), status = c(status_vec, 0))
            }
          } else {
            list(residual = final_val, status = 0)
          }
        },
        by = .(pair_id)
      ]
    }
  }

  resid_data_1_0 <- NULL
  if (get_1_0 && !is.null(model_1_0)) {
    data_1_0 <- if (!is.null(model_1_0$data)) {
      model_1_0$data
    } else if (!is.null(reproduced_data)) {
      reproduced_data$data_1_0
    } else {
      NULL
    }

    if (!is.null(data_1_0)) {
      dt_1_0 <- data.table::as.data.table(data_1_0)

      resid_data_1_0 <- dt_1_0[,
        {
          # Access the prediction column if exists, otherwise compute it
          curr_pred_vals <- if ("prediction" %in% names(.SD)) .SD$prediction else calculate_predictions_helper(model_1_0, .SD)

          cum_lambda <- cumsum(as.double(curr_pred_vals))
          event_indices <- which(event == 1)
          final_val <- if (length(cum_lambda) > 0) tail(cum_lambda, 1) else 0

          if (length(event_indices) > 0) {
            event_vals <- cum_lambda[event_indices]
            event_resids <- diff(c(0, event_vals))
            status_vec <- rep(1, length(event_resids))

            tail_resid <- final_val - tail(event_vals, 1)
            if (abs(tail_resid) < 1e-10) {
              list(residual = event_resids, status = status_vec)
            } else {
              list(residual = c(event_resids, tail_resid), status = c(status_vec, 0))
            }
          } else {
            list(residual = final_val, status = 0)
          }
        },
        by = .(pair_id)
      ]
    }
  }

  raw_list <- list(resid_0_1 = resid_data_0_1, resid_1_0 = resid_data_1_0)

  if (raw) {
    class(raw_list) <- "redeem_residuals"
    return(raw_list)
  }

  # Default: Return KM survival estimates frame
  return(prepare_residual_plot_data(raw_list))
}

#' Internal helper to prepare residual plot data
#' @param raw_list List of raw residuals
#' @importFrom survival Surv survfit
#' @keywords internal
prepare_residual_plot_data <- function(raw_list) {
  # Filter out NULL entries (e.g., resid_1_0 for REM or when get_* is FALSE)
  valid_list <- Filter(Negate(is.null), raw_list)
  if (length(valid_list) == 0) {
    return(NULL)
  }

  combined <- data.table::rbindlist(valid_list)
  if (nrow(combined) == 0) {
    return(NULL)
  }

  # Suppress warnings from survfit about 0 time if any (though unlikely here)
  fit <- survival::survfit(survival::Surv(combined$residual, combined$status) ~ 1)

  res <- data.frame(
    time = fit$time,
    surv = fit$surv,
    lower = fit$lower,
    upper = fit$upper,
    theoretical = exp(-fit$time)
  )
  class(res) <- c("redeem_residual_plot_data", class(res))
  return(res)
}
