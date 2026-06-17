#' @importFrom graphics par plot lines title mtext
#' @importFrom grDevices rainbow
NULL

#' Print DEM Model
#'
#' @param x A `dem` object.
#' @param ... Additional arguments.
#' @method print dem
#' @noRd
#' @export
print.dem <- function(x, ...) {
  cat("Durational Event Model (DEM)\n")
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Number of nodes:", x$n_nodes, "\n")
  cat("Number of events:", sum(x$event_numbers), "\n")
  if (!is.null(x$runtime)) {
    cat("Estimation time:", format(x$runtime), "\n")
  }
  invisible(x)
}

#' Summary of DEM Model Results
#'
#' @param object A `dem` object.
#' @param ... Additional arguments.
#' @method summary dem
#' @noRd
#' @export
summary.dem <- function(object, ...) {
  aic_0_1 <- aic_1_0 <- bic_0_1 <- bic_1_0 <- NA

  tmp_est_incidence <- NULL
  if (!is.null(object$model_0_1)) {
    tmp_est_incidence <- summary(object$model_0_1)
    if(inherits(object$model_0_1, "coxph")){
      tmp_est_incidence$coefficients <- round(tmp_est_incidence$coefficients, digits = 3)
    }
  }

  tmp_est_duration <- NULL
  if (!is.null(object$model_1_0)) {
    tmp_est_duration <- summary(object$model_1_0)
    if(inherits(object$model_1_0, "coxph")){
      tmp_est_duration$coefficients <- round(tmp_est_duration$coefficients, digits = 3)
    }
  }

  if (!is.null(object$model_0_1) && !is.null(object$model_1_0)) {
    if (inherits(object$model_1_0, "dem.mm")) {
      aic_0_1 <- 2 * (length(object$model_0_1$est_core) + length(object$model_0_1$est_time)+ length(object$model_0_1$est_degree))- 2 * object$model_0_1$llh
      aic_1_0 <- 2 * (length(object$model_1_0$est_core) + length(object$model_1_0$est_time)+ length(object$model_1_0$est_degree))- 2 * object$model_1_0$llh
      bic_0_1 <- 2 * log(object$event_numbers[1]) - 2 * object$model_0_1$llh
      bic_1_0 <- 2 * log(object$event_numbers[2]) - 2 * object$model_1_0$llh
    } else {
      aic_0_1 <- 2 * (length(object$model_0_1$est_core) + length(object$model_0_1$est_time)+ length(object$model_0_1$est_degree)) - 2 * logLik(object$model_0_1)[1]
      aic_1_0 <- 2 * (length(object$model_1_0$est_core) + length(object$model_1_0$est_time)+ length(object$model_1_0$est_degree)) - 2 * logLik(object$model_1_0)[1]
      bic_0_1 <- 2 * log(object$event_numbers[1]) - 2 * logLik(object$model_0_1)[1]
      bic_1_0 <- 2 * log(object$event_numbers[2]) - 2 * logLik(object$model_1_0)[1]
    }
    AIC <- aic_0_1 + aic_1_0
    BIC <- bic_0_1 + bic_1_0
  } else {
    AIC <- NA
    BIC <- NA
  }

  res <- list(
    call = object$call,
    Incidence = tmp_est_incidence,
    Duration = tmp_est_duration,
    AIC = AIC,
    aic_0_1 = aic_0_1,
    aic_1_0 = aic_1_0,
    BIC = BIC,
    runtime = object$runtime
  )
  class(res) <- "summary.dem"
  return(res)
}

#' Print Summary of DEM Results
#'
#' @param x A `summary.dem` object.
#' @param ... Additional arguments.
#' @method print summary.dem
#' @noRd
#' @export
print.summary.dem <- function(x, ...) {
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  if (!is.null(x$Incidence)) {
    cat("Results for Incidence Intensity (0 -> 1): \n")
    print(x$Incidence, ...)
  }

  if (!is.null(x$Duration)) {
    cat("\nResults for Duration Intensity (1 -> 0): \n")
    print(x$Duration, ...)
  }

  if (!is.na(x$AIC)) {
    cat("\nCombined Model Fit:\n")
    cat("  AIC:", x$AIC, "\n  BIC:", x$BIC, "\n")
  }

  if (!is.null(x$runtime)) {
    cat("\nTotal estimation time:", format(x$runtime), "\n")
  }
  invisible(x)
}

#' Summary of a \code{redeem_result} Model Fit
#'
#' Computes a summary of a fitted \code{redeem_result} object, collecting
#' the estimated fixed effects, log-likelihood, and (if present) degree and
#' temporal baseline effects into a structured list suitable for printing.
#'
#' @method summary redeem_result
#' @param object A `redeem_result` object.
#' @param ... Additional arguments (currently unused).
#' @return An object of class \code{summary.redeem_result}, which is a list
#'   containing:
#' \itemize{
#'   \item \code{coefficients}: A numeric matrix with one row per fixed-effect
#'     covariate and columns \code{Estimate}, \code{Std. Error}, \code{t value},
#'     and \code{Pr(>|t|)}.
#'   \item \code{llh}: The log-likelihood of the fitted model (scalar).
#'   \item \code{degree_summary}: A list with summary statistics (\code{n},
#'     \code{n_unidentifiable}, \code{mean}, \code{sd}, \code{range}) of the
#'     estimated degree effects, only present when more than 10 degree
#'     parameters were estimated.
#'   \item \code{degree_effects}: A named numeric vector of estimated degree
#'     effects, only present when 10 or fewer degree parameters were estimated.
#'   \item \code{time_summary}: A list with summary statistics of the estimated
#'     temporal baseline effects, only present when more than 10 time intervals
#'     were used.
#'   \item \code{time_effects}: A named numeric vector of estimated temporal
#'     baseline effects, only present when 10 or fewer time intervals were used.
#'   \item \code{iter}: Integer; the number of iterations performed by the
#'     optimizer (\code{NA} if history was not saved).
#' }
#' @export
summary.redeem_result <- function(object, ...) {
  # Core effects
  est <- object$est_core

  # Calculate Std. Errors for core effects if available
  stderr <- rep(NA, length(est))
  if (!is.null(object$covariance) && length(est) > 0) {
    diag_vals <- diag(object$covariance)
    diag_vals[diag_vals < 0] <- 0
    stderr[seq_along(est)] <- sqrt(diag_vals)
  }

  tvalue <- est / stderr
  # Handle cases where stderr is 0 or NA
  tvalue[is.nan(tvalue)] <- NA
  pvalue <- 2 * stats::pnorm(-abs(tvalue))
  coef_table <- cbind(est, stderr, tvalue, pvalue)
  colnames(coef_table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  if (!is.null(names(est))) rownames(coef_table) <- names(est)

  res <- list(coefficients = coef_table, llh = object$llh)

  # Degree Effects Summary
  if (!is.null(object$est_degree) && length(object$est_degree) > 0) {
    unidentifiable <- is.infinite(object$est_degree) & (object$est_degree < 0)
    finite_degrees <- object$est_degree[!unidentifiable]
    if (length(object$est_degree) > 10) {
      res$degree_summary <- list(
        n = length(object$est_degree),
        n_unidentifiable = sum(unidentifiable),
        mean = mean(finite_degrees, na.rm = TRUE),
        sd = stats::sd(finite_degrees, na.rm = TRUE),
        range = if (length(finite_degrees) > 0) range(finite_degrees, na.rm = TRUE) else c(NA, NA)
      )
    } else {
      res$degree_effects <- object$est_degree
    }
  }

  # Temporal Effects Summary
  if (!is.null(object$est_time) && length(object$est_time) > 0) {
    # Check if they are already in the coefficients table names
    in_coefs <- all(names(object$est_time) %in% rownames(coef_table))

    unidentifiable_t <- is.infinite(object$est_time) & (object$est_time < 0)
    finite_time <- object$est_time[!unidentifiable_t]

    if (length(object$est_time) > 10) {
      res$time_summary <- list(
        n = length(object$est_time),
        n_unidentifiable = sum(unidentifiable_t),
        mean = mean(finite_time, na.rm = TRUE),
        sd = stats::sd(finite_time, na.rm = TRUE),
        range = if (length(finite_time) > 0) range(finite_time, na.rm = TRUE) else c(NA, NA)
      )
    } else if (!in_coefs) {
      res$time_effects <- object$est_time
    }
  }

  res$iter <- if (!is.null(object$llh_hist)) length(object$llh_hist) else NA

  class(res) <- "summary.redeem_result"
  return(res)
}

#' @method print summary.redeem_result
#' @noRd
#' @export
print.summary.redeem_result <- function(x, ...) {
  if (!is.null(x$coefficients) && nrow(x$coefficients) > 0) {
    cat("Fixed Effects:\n")
    stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE, ...)
  } else {
    cat("No fixed covariate effects.\n")
  }

  if (!is.null(x$degree_summary)) {
    cat("\nDegree Effects Summary:\n")
    cat("  Nodes:", x$degree_summary$n, " (Unidentifiable:", x$degree_summary$n_unidentifiable, ")\n")
    cat("  Mean:", round(x$degree_summary$mean, 4), " SD:", round(x$degree_summary$sd, 4), "\n")
    cat("  Range: [", round(x$degree_summary$range[1], 4), ", ", round(x$degree_summary$range[2], 4), "]\n")
  } else if (!is.null(x$degree_effects)) {
    cat("\nDegree Effects:\n")
    print(round(x$degree_effects, 4))
  }

  if (!is.null(x$time_summary)) {
    cat("\nTemporal Effects Summary:\n")
    cat("  Intervals:", x$time_summary$n, " (Unidentifiable:", x$time_summary$n_unidentifiable, ")\n")
    cat("  Mean:", round(x$time_summary$mean, 4), " SD:", round(x$time_summary$sd, 4), "\n")
    cat("  Range: [", round(x$time_summary$range[1], 4), ", ", round(x$time_summary$range[2], 4), "]\n")
  } else if (!is.null(x$time_effects)) {
    cat("\nTemporal Effects:\n")
    print(round(x$time_effects, 4))
  }

  if (!is.null(x$llh)) {
    cat("\nLog-likelihood:", round(x$llh, 3), "\n")
  }
  invisible(x)
}


#' Log-likelihood for redeem_result object
#'
#' @param object A `redeem_result` object.
#' @param ... Additional arguments.
#' @method logLik redeem_result
#' @noRd
#' @export
logLik.redeem_result <- function(object, ...) {
  res <- object$llh
  n_finite <- sum(is.finite(object$est_core)) + sum(is.finite(object$est_degree)) + sum(is.finite(object$est_time))
  attr(res, "df") <- n_finite
  attr(res, "nobs") <- object$n_obs
  class(res) <- "logLik"
  return(res)
}

#' Print REM Model
#'
#' @param x A `rem` object.
#' @param ... Additional arguments.
#' @method print rem
#' @noRd
#' @export
print.rem <- function(x, ...) {
  cat("Relational Event Model (REM)\n")
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Number of nodes:", x$n_nodes, "\n")
  cat("Number of events:", sum(x$event_numbers), "\n")
  if (!is.null(x$runtime)) {
    cat("Estimation time:", format(x$runtime), "\n")
  }
  invisible(x)
}

#' Summary of REM Model Results
#'
#' @param object A `rem` object.
#' @param ... Additional arguments.
#' @method summary rem
#' @noRd
#' @export
summary.rem <- function(object, ...) {
  res <- summary(object$model)
  res$call <- object$call
  res$runtime <- object$runtime
  class(res) <- c("summary.rem", class(res))
  return(res)
}

#' Print Summary of REM Results
#'
#' @param x A `summary.rem` object.
#' @param ... Additional arguments.
#' @method print summary.rem
#' @noRd
#' @export
print.summary.rem <- function(x, ...) {
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  # Exclude summary.rem from class to delegate to print.summary.redeem_result
  class(x) <- setdiff(class(x), "summary.rem")
  print(x, ...)
  if (!is.null(x$runtime)) {
    cat("\nEstimation time:", format(x$runtime), "\n")
  }
  invisible(x)
}

#' Print Ranking Results
#'
#' @param x A `ranking_redeem` object.
#' @param ... Additional arguments.
#' @method print ranking_redeem
#' @noRd
#' @export
print.ranking_redeem <- function(x, ...) {
  cat("Ranking Results (Redeem)\n")
  mrr <- attr(x, "mrr")
  if (!is.null(mrr)) {
    cat("Mean Reciprocal Rank (MRR):", mrr, "\n")
  }
  mean_rank <- attr(x, "mean_rank")
  if (!is.null(mean_rank)) {
    cat("Mean Rank:", mean_rank, "\n")
  }
  median_rank <- attr(x, "median_rank")
  if (!is.null(median_rank)) {
    cat("Median Rank:", median_rank, "\n")
  }
  hits_summary <- attr(x, "hits_summary")
  if (!is.null(hits_summary)) {
    cat("\nTop-K Goodness-of-Fit Summary:\n")
    print(hits_summary)
    cat("\n")
  }
  print(head(as.data.frame(x), 10))
  if (nrow(x) > 10) cat("... (truncated)\n")
  invisible(x)
}

#' Log-likelihood for rem object
#'
#' @param object A `rem` object.
#' @param ... Additional arguments.
#' @method logLik rem
#' @noRd
#' @export
logLik.rem <- function(object, ...) {
  return(logLik(object$model))
}

#' Plot Convergence of redeem_result model
#'
#' @param x A `redeem_result` object.
#' @param coefs Logical; plot fixed coefficients.
#' @param degree Logical; plot degree effects.
#' @param time Logical; plot temporal effects.
#' @param llh Logical; plot log-likelihood history.
#' @param sub_label Character; subtitle to display at the bottom of the figure.
#' @param separate Logical; if `TRUE`, each plot is generated in a new window
#'   or as a separate sequence.
#' @param baseline Logical; if `TRUE`, plot the estimated baseline intensity.
#' @param ... Additional arguments passed to `graphics::plot`.
#' @method plot redeem_result
#' @noRd
#' @export
plot.redeem_result <- function(x, ..., coefs = TRUE, degree = TRUE,
                             time = TRUE, llh = TRUE, baseline = FALSE,
                             sub_label = NULL, separate = FALSE) {
  # Save original par settings that we might change and ensure they are restored
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  # Build plot indicators based on split histories if available
  plot_fixed <- coefs && !is.null(x$coefficients_core_hist) && ncol(x$coefficients_core_hist) > 0
  plot_degree <- degree && !is.null(x$est_degree_hist) && ncol(x$est_degree_hist) > 0
  plot_time <- time && !is.null(x$est_time_hist) && ncol(x$est_time_hist) > 0

  # Fallback to combined coef_hist if split histories are missing (backward compatibility)
  if (!plot_fixed && !plot_degree && !plot_time && !is.null(x$coef_hist)) {
    nm <- colnames(x$coef_hist)
    is_degree <- grep("^effect_|^sender_|^receiver_", nm)
    is_time <- grep("^time_cat", nm)
    is_fixed <- setdiff(seq_along(nm), c(is_degree, is_time))

    plot_fixed <- coefs && length(is_fixed) > 0
    plot_degree <- degree && length(is_degree) > 0
    plot_time <- time && length(is_time) > 0
  } else {
    is_fixed <- is_degree <- is_time <- NULL
  }

  # Check if we were expecting plots but history was suppressed
  if ((coefs && !plot_fixed) || (degree && !plot_degree) || (time && !plot_time)) {
    if (is.null(x$coefficients_core_hist) && is.null(x$est_degree_hist) && is.null(x$coef_hist)) {
        message("Trace plots for coefficients are unavailable because history was not saved during estimation.")
        message("To enable these plots, set 'save_hist = TRUE' in control.redeem().")
    }
  }

  n_plots <- sum(c(plot_fixed, plot_degree, plot_time, llh && !is.null(x$llh_hist), baseline))
  if (n_plots == 0) {
    if (!llh || is.null(x$llh_hist)) {
        message("No history data available to plot.")
    }
    return(invisible(x))
  }
  if (!separate) {
    if (n_plots > 1 || !is.null(sub_label)) {
      # Reserve space at the bottom for the sub_label
      oma <- if (!is.null(sub_label)) c(4, 0, 0, 0) else c(0, 0, 0, 0)
      graphics::par(mfrow = c(ceiling(n_plots / 2), 2), oma = oma, las = 1, bty = "l")
    } else {
      graphics::par(las = 1, bty = "l")
    }
  } else {
    graphics::par(las = 1, bty = "l")
  }

  if (llh && !is.null(x$llh_hist)) {
    main_llh <- "Log-likelihood Trace"
    if (separate && !is.null(sub_label)) main_llh <- paste0(main_llh, " (", sub_label, ")")
    graphics::plot(x$llh_hist, type = "l", xlab = "Iteration", ylab = "Log-likelihood",
         main = main_llh)
  }

  if (plot_fixed) {
    hist_f <- if (!is.null(is_fixed)) x$coef_hist[, is_fixed, drop = FALSE] else x$coefficients_core_hist
    main_fixed <- "Core Effects (Relational)"
    if (separate && !is.null(sub_label)) main_fixed <- paste0(main_fixed, " (", sub_label, ")")
    graphics::plot(NA, xlim = c(1, nrow(hist_f)),
         ylim = range(hist_f, na.rm = TRUE),
         xlab = "Iteration", ylab = "Estimates",
         main = main_fixed)
    cols <- grDevices::rainbow(ncol(hist_f))
    for (i in seq_len(ncol(hist_f))) {
      graphics::lines(y = hist_f[, i], x = seq_len(nrow(hist_f)), col = cols[i])
    }
  }

  if (plot_degree) {
    hist_p <- if (!is.null(is_degree)) x$coef_hist[, is_degree, drop = FALSE] else x$est_degree_hist
    main_degree <- "Degree Effects"
    if (separate && !is.null(sub_label)) main_degree <- paste0(main_degree, " (", sub_label, ")")
    graphics::plot(NA, xlim = c(1, nrow(hist_p)),
         ylim = range(hist_p, na.rm = TRUE),
         xlab = "Iteration", ylab = "Estimates",
         main = main_degree)
    cols <- grDevices::rainbow(ncol(hist_p))
    for (i in seq_len(ncol(hist_p))) {
      graphics::lines(y = hist_p[, i], x = seq_len(nrow(hist_p)), col = cols[i])
    }
  }

  if (plot_time) {
    hist_t <- if (!is.null(is_time)) x$coef_hist[, is_time, drop = FALSE] else x$est_time_hist
    main_time <- "Temporal Effects"
    if (separate && !is.null(sub_label)) main_time <- paste0(main_time, " (", sub_label, ")")
    graphics::plot(NA, xlim = c(1, nrow(hist_t)),
         ylim = range(hist_t, na.rm = TRUE),
         xlab = "Iteration", ylab = "Estimates",
         main = main_time)
    cols <- grDevices::rainbow(ncol(hist_t))
    for (i in seq_len(ncol(hist_t))) {
      graphics::lines(y = hist_t[, i], x = seq_len(nrow(hist_t)), col = cols[i])
    }
  }
  if (baseline) {
    plot_baseline(x, ..., separate = separate, sub_label = sub_label)
  }

  if (!separate && !is.null(sub_label)) {
    graphics::mtext(sub_label, side = 1, outer = TRUE, line = 1.5, font = 2, cex = 1.2)
  }
  invisible(x)
}

#' Plot Trace Plots for DEM Model
#'
#' @param x A `dem` object.
#' @param which Integer indicating which transition plots to show (1: 0->1,
#'   2: 1->0, 3: both).
#' @param separate Logical; plot each component separately.
#' @param baseline Logical; plot the estimated baseline intensity.
#' @param ... Additional arguments passed to specific plot methods.
#' @method plot dem
#' @noRd
#' @export
plot.dem <- function(x, which = 3, separate = FALSE, baseline = FALSE, ...) {
  # No local par changes needed here as plot.redeem_result handles its own state

  if (which %in% c(1, 3) && !is.null(x$model_0_1)) {
    plot(x$model_0_1, ..., separate = separate, baseline = baseline, sub_label = "Transition 0 -> 1")
  }
  if (which %in% c(2, 3) && !is.null(x$model_1_0)) {
    plot(x$model_1_0, ..., separate = separate, baseline = baseline, sub_label = "Transition 1 -> 0")
  }
}

#' Plot Trace Plots for REM Model
#'
#' @param x A `rem` object.
#' @param baseline Logical; plot the estimated baseline intensity.
#' @param ... Additional arguments passed to specific plot methods.
#' @method plot rem
#' @noRd
#' @export
plot.rem <- function(x, baseline = FALSE, ...) {
  if (!is.null(x$model)) {
    plot(x$model, ..., baseline = baseline)
  } else {
    message("No model found in rem object.")
  }
}

#' Plot the Estimated Baseline Intensity
#'
#' Draws a step-function plot of the estimated piecewise-constant baseline
#' intensity against time. The function dispatches to class-specific methods
#' for \code{\link{dem}}, \code{\link{rem}}, and \code{redeem_result} objects.
#'
#' @param x A \code{\link{dem}}, \code{\link{rem}}, or
#'   \code{redeem_result} object produced by \code{\link{dem}} or
#'   \code{\link{rem}}.
#' @param ... Additional arguments passed to \code{graphics::plot}.
#' @return The original object \code{x} is returned invisibly. Called
#'   primarily for its side effect of producing a plot.
#' @export
plot_baseline <- function(x, ...) {
  UseMethod("plot_baseline")
}

#' Plot Baseline for DEM
#'
#' @param x A `dem` object.
#' @method plot_baseline dem
#' @param process Character; either "formation" (0->1) or "dissolution" (1->0).
#' @param ... Additional arguments passed to [graphics::plot()].
#' @noRd
#' @export
plot_baseline.dem <- function(x, process = c("formation", "dissolution"), ...) {
  process <- match.arg(process)
  if (process == "formation") {
    if (is.null(x$model_0_1)) stop("Formation model was not estimated.")
    plot_baseline(x$model_0_1, main = "Baseline Intensity: Formation", ...)
  } else {
    if (is.null(x$model_1_0)) stop("Dissolution model was not estimated.")
    plot_baseline(x$model_1_0, main = "Baseline Intensity: Dissolution", ...)
  }
}

#' Plot Baseline for REM
#'
#' @param x A `rem` object.
#' @method plot_baseline rem
#' @param ... Additional arguments passed to [graphics::plot()].
#' @noRd
#' @export
plot_baseline.rem <- function(x, ...) {
  if (is.null(x$model)) stop("REM model was not estimated.")
  plot_baseline(x$model, main = "Baseline Intensity", ...)
}

#' Plot Baseline for redeem_result
#'
#' @param x A `redeem_result` object.
#' @method plot_baseline redeem_result
#' @param sub_label Character; subtitle for the plot.
#' @param separate Logical; if `TRUE`, uses `par` settings for a standalone plot.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @noRd
#' @export
plot_baseline.redeem_result <- function(x, ..., sub_label = NULL, separate = TRUE) {
  if (is.null(x$est_time)) {
    message("No temporal baseline effects found in the model.")
    return(invisible(x))
  }

  if (separate) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(las = 1, bty = "l")
  }

  changepoints <- x$time_changepoints
  est_time <- x$est_time
  full_baseline <- x$full_baseline

  # Intensities
  intensities <- if (full_baseline) exp(est_time) else exp(c(0, est_time))

  # Time intervals
  # The intervals are [0, t1), [t1, t2), ..., [tM, max_time]
  max_time <- if (!is.null(x$data)) {
    m_val <- max(x$data$time_new, na.rm = TRUE)
    if (!is.finite(m_val)) m_val <- max(changepoints, 0, na.rm = TRUE) + 1
    m_val
  } else {
    max(changepoints, 0, na.rm = TRUE) + 1
  }

  times <- c(0, changepoints, max_time)

  # Prepare title
  main_title <- if (!is.null(list(...)$main)) list(...)$main else "Baseline Intensity"
  if (separate && !is.null(sub_label)) main_title <- paste0(main_title, " (", sub_label, ")")

  # Plot step function
  plot_args <- list(x = times, y = c(intensities, tail(intensities, 1)),
                    type = "s", main = main_title, xlab = "Time", ylab = "Intensity")
  user_args <- list(...)
  # User arguments override defaults
  for (n in names(user_args)) plot_args[[n]] <- user_args[[n]]

  do.call(graphics::plot, plot_args)

  invisible(x)
}

#' Plot Residuals Diagnostic
#'
#' @param x A `redeem_residual_plot_data` or `redeem_residuals` object.
#' @param ... Additional arguments passed to `graphics::plot`.
#' @method plot redeem_residual_plot_data
#' @noRd
#' @export
plot.redeem_residual_plot_data <- function(x, ...) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(las = 1, bty = "l")

  graphics::plot(x$time, x$surv, type = "l",
                 xlab = "Residuals", ylab = "Survival Probability", ...)
  graphics::lines(x$time, x$theoretical, col = "#34A853", lty = 2)
  graphics::legend("topright", legend = c("Estimated", "Theoretical (Exp(1))"),
                   col = c("black", "#34A853"), lty = c(1, 2), bty = "n")
  invisible(x)
}

#' Plot Raw Residuals
#'
#' @param x A `redeem_residuals` object.
#' @param ... Additional arguments passed to `graphics::plot`.
#' @method plot redeem_residuals
#' @noRd
#' @export
plot.redeem_residuals <- function(x, ...) {
  # Convert raw residuals to plot data
  plot_data <- prepare_residual_plot_data(x)
  if (is.null(plot_data)) {
    message("No valid residuals to plot.")
    return(invisible(x))
  }
  plot(plot_data, ...)
}

#' Plot Ranking Goodness-of-Fit
#'
#' @param x A `ranking_redeem` object.
#' @param metric Character indicating which metric to plot: "recall",
#'   "precision", or "both".
#' @param ... Additional arguments passed to `graphics::plot`.
#' @method plot ranking_redeem
#' @noRd
#' @export
plot.ranking_redeem <- function(x, metric = c("recall", "precision", "both"), ...) {
  metric <- match.arg(metric)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(las = 1, bty = "l")

  if (metric == "recall") {
    graphics::plot(x$Cutpoint, x$Recall, type = "l", col = "blue", lwd = 2,
                   xlab = "Cutpoint (k)", ylab = "Recall", main = "Recall GOF", ...)
  } else if (metric == "precision") {
    graphics::plot(x$Cutpoint, x$Precision, type = "l", col = "red", lwd = 2,
                   xlab = "Cutpoint (k)", ylab = "Precision", main = "Precision GOF", ...)
  } else {
    # metric == "both"
    graphics::plot(x$Cutpoint, x$Recall, type = "l", col = "blue", lwd = 2,
                   xlab = "Cutpoint (k)", ylab = "Metric Value", main = "GOF Metrics",
                   ylim = c(0, max(c(x$Recall, x$Precision), na.rm = TRUE)), ...)
    graphics::lines(x$Cutpoint, x$Precision, col = "red", lwd = 2, lty = 2)
    graphics::legend("topleft", legend = c("Recall", "Precision"),
                     col = c("blue", "red"), lty = c(1, 2), lwd = 2)
  }
  invisible(x)
}


#' @method predict rem
#' @noRd
#' @export
predict.rem <- function(object, time = NULL, type = c("response", "lp", "terms"), ...) {
  type <- match.arg(type)

  # Allow for the function to use the data object if available, else generate it on the spot
  data <- if (!is.null(object$model$data)) {
    object$model$data
  } else {
    reproduce_model_data(object)
  }

  if (!is.null(time)) {
    # Filter to intervals containing any of the specified time points
    keep <- rowSums(sapply(time, function(t) data$time <= t & data$time_new > t)) > 0
    data <- data[keep, ]
  }

  res <- predict_transition_helper(object$model, data, type)

  from_vec <- if (!is.null(data) && !is.null(data$from)) data$from else integer(0)
  to_vec <- if (!is.null(data) && !is.null(data$to)) data$to else integer(0)

  if (type %in% c("response", "lp")) {
    df <- data.frame(
      from = from_vec,
      to = to_vec,
      prediction = as.numeric(res),
      mode = rep("formation", length(res)),
      stringsAsFactors = FALSE
    )
    return(df)
  } else {
    df <- data.frame(
      from = from_vec,
      to = to_vec,
      mode = rep("formation", nrow(res)),
      as.data.frame(res),
      stringsAsFactors = FALSE
    )
    return(df)
  }
}

#' @method predict dem
#' @noRd
#' @export
predict.dem <- function(object, time = NULL, type = c("response", "lp", "terms"), process = c("both", "formation", "dissolution"), ...) {
  type <- match.arg(type)
  process <- match.arg(process)

  # Allow for the function to use the data object if available, else generate it on the spot
  data_list <- if (!is.null(object$model_0_1$data) || !is.null(object$model_1_0$data)) {
    list(
      data_0_1 = object$model_0_1$data,
      data_1_0 = object$model_1_0$data
    )
  } else {
    reproduce_model_data(object)
  }

  # 1. Formation (0 -> 1)
  df_form <- NULL
  if (process %in% c("both", "formation") && !is.null(object$model_0_1)) {
    data_0_1 <- data_list$data_0_1
    if (!is.null(data_0_1) && nrow(data_0_1) > 0) {
      if (!is.null(time)) {
        keep <- rowSums(sapply(time, function(t) data_0_1$time <= t & data_0_1$time_new > t)) > 0
        data_0_1 <- data_0_1[keep, ]
      }
    }
    if (!is.null(data_0_1) && nrow(data_0_1) > 0) {
      pred_form <- predict_transition_helper(object$model_0_1, data_0_1, type)
      from_vec <- if (!is.null(data_0_1$from)) data_0_1$from else integer(0)
      to_vec <- if (!is.null(data_0_1$to)) data_0_1$to else integer(0)
      if (type %in% c("response", "lp")) {
        df_form <- data.frame(
          from = from_vec,
          to = to_vec,
          prediction = as.numeric(pred_form),
          mode = "formation",
          stringsAsFactors = FALSE
        )
      } else {
        df_form <- data.frame(
          from = from_vec,
          to = to_vec,
          mode = "formation",
          as.data.frame(pred_form),
          stringsAsFactors = FALSE
        )
      }
    } else {
      pred_form <- predict_transition_helper(object$model_0_1, NULL, type)
      if (type %in% c("response", "lp")) {
        df_form <- data.frame(
          from = integer(0),
          to = integer(0),
          prediction = numeric(0),
          mode = character(0),
          stringsAsFactors = FALSE
        )
      } else {
        df_form <- data.frame(
          from = integer(0),
          to = integer(0),
          mode = character(0),
          as.data.frame(pred_form),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # 2. Dissolution (1 -> 0)
  df_diss <- NULL
  if (process %in% c("both", "dissolution") && !is.null(object$model_1_0)) {
    data_1_0 <- data_list$data_1_0
    if (!is.null(data_1_0) && nrow(data_1_0) > 0) {
      if (!is.null(time)) {
        keep <- rowSums(sapply(time, function(t) data_1_0$time <= t & data_1_0$time_new > t)) > 0
        data_1_0 <- data_1_0[keep, ]
      }
    }
    if (!is.null(data_1_0) && nrow(data_1_0) > 0) {
      pred_diss <- predict_transition_helper(object$model_1_0, data_1_0, type)
      from_vec <- if (!is.null(data_1_0$from)) data_1_0$from else integer(0)
      to_vec <- if (!is.null(data_1_0$to)) data_1_0$to else integer(0)
      if (type %in% c("response", "lp")) {
        df_diss <- data.frame(
          from = from_vec,
          to = to_vec,
          prediction = as.numeric(pred_diss),
          mode = "dissolution",
          stringsAsFactors = FALSE
        )
      } else {
        df_diss <- data.frame(
          from = from_vec,
          to = to_vec,
          mode = "dissolution",
          as.data.frame(pred_diss),
          stringsAsFactors = FALSE
        )
      }
    } else {
      pred_diss <- predict_transition_helper(object$model_1_0, NULL, type)
      if (type %in% c("response", "lp")) {
        df_diss <- data.frame(
          from = integer(0),
          to = integer(0),
          prediction = numeric(0),
          mode = character(0),
          stringsAsFactors = FALSE
        )
      } else {
        df_diss <- data.frame(
          from = integer(0),
          to = integer(0),
          mode = character(0),
          as.data.frame(pred_diss),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (process == "formation") {
    return(df_form)
  } else if (process == "dissolution") {
    return(df_diss)
  } else {
    # process == "both"
    if (type %in% c("response", "lp")) {
      combined_df <- rbind(df_form, df_diss)
      return(combined_df)
    } else {
      # type == "terms"
      if (!is.null(df_form) && ncol(df_form) > 3) {
        term_cols <- setdiff(names(df_form), c("from", "to", "mode"))
        names(df_form)[names(df_form) %in% term_cols] <- paste0("formation_", term_cols)
      }
      if (!is.null(df_diss) && ncol(df_diss) > 3) {
        term_cols <- setdiff(names(df_diss), c("from", "to", "mode"))
        names(df_diss)[names(df_diss) %in% term_cols] <- paste0("dissolution_", term_cols)
      }

      combined_dt <- data.table::rbindlist(list(df_form, df_diss), fill = TRUE)

      term_cols_combined <- setdiff(names(combined_dt), c("from", "to", "mode"))
      for (col in term_cols_combined) {
        data.table::set(combined_dt, which(is.na(combined_dt[[col]])), col, 0)
      }

      return(as.data.frame(combined_dt))
    }
  }
}

#' Helper to predict single transition process
#' @keywords internal
predict_transition_helper <- function(model, data, type) {
  nm <- names(model$coefficients)
  is_degree <- grep("^effect_|^sender_|^receiver_", nm)
  is_time <- grep("time|duration|elapsed|baseline", nm, ignore.case = TRUE)
  is_core <- setdiff(seq_along(nm), c(is_degree, is_time))

  if (is.null(data) || nrow(data) == 0) {
    if (type == "terms") {
      # Build an empty matrix with the structurally correct column names
      term_names <- c()
      if (length(is_core) > 0) term_names <- c(term_names, nm[is_core])
      if (!is.null(model$est_degree) && length(model$est_degree) > 0) {
        directed <- model$directed
        if (is.null(directed)) {
          directed <- length(model$est_degree) > 0 # Default/fallback
        }
        if (isTRUE(directed)) {
          term_names <- c(term_names, "sender_effects", "receiver_effects")
        } else {
          term_names <- c(term_names, "degree_effects")
        }
      }
      if ((!is.null(model$est_time) && length(model$est_time) > 0 && !is.null(model$time_changepoints)) || length(is_time) > 0) {
        term_names <- c(term_names, "time_effects")
      }
      empty_mat <- matrix(0, nrow = 0, ncol = length(term_names))
      colnames(empty_mat) <- term_names
      return(empty_mat)
    } else {
      return(numeric(0))
    }
  }

  # 1. Compute Core terms contributions
  core_contribs <- list()
  lp_core <- rep(0, nrow(data))
  if (length(is_core) > 0) {
    for (c_name in nm[is_core]) {
      contrib <- rep(0, nrow(data))
      if (c_name %in% names(data)) {
        contrib <- data[[c_name]] * model$coefficients[c_name]
      } else if (tolower(c_name) %in% c("intercept", "(intercept)")) {
        # Case-insensitive robust fallback for intercept
        matching_names <- names(data)[tolower(names(data)) == tolower(c_name)]
        if (length(matching_names) > 0) {
          contrib <- data[[matching_names[1]]] * model$coefficients[c_name]
        } else {
          # Default to intercept constant of 1
          contrib <- rep(1, nrow(data)) * model$coefficients[c_name]
        }
      }
      core_contribs[[c_name]] <- contrib
      lp_core <- lp_core + contrib
    }
  }

  # 2. Compute Degree effects contributions
  deg_contribs <- list()
  lp_degree <- rep(0, nrow(data))
  deg_vals <- model$est_degree
  if (!is.null(deg_vals) && length(deg_vals) > 0) {
    directed <- model$directed
    if (is.null(directed)) {
      observed_n_nodes <- max(c(data$from, data$to), na.rm = TRUE)
      directed <- length(deg_vals) >= 2 * observed_n_nodes
    }
    n_nodes <- model$n_nodes
    if (is.null(n_nodes)) {
      observed_n_nodes <- max(c(data$from, data$to), na.rm = TRUE)
      if (directed) {
        n_nodes <- min(observed_n_nodes, floor(length(deg_vals) / 2))
      } else {
        n_nodes <- min(observed_n_nodes, length(deg_vals))
      }
    }
    if (directed) {
      mu <- deg_vals[1:n_nodes]
      nu <- deg_vals[(n_nodes + 1):(2 * n_nodes)]
      sender_contrib <- mu[data$from]
      receiver_contrib <- nu[data$to]

      sender_contrib[is.na(sender_contrib)] <- 0
      receiver_contrib[is.na(receiver_contrib)] <- 0

      deg_contribs[["sender_effects"]] <- sender_contrib
      deg_contribs[["receiver_effects"]] <- receiver_contrib
      lp_degree <- sender_contrib + receiver_contrib
    } else {
      degree_contrib <- deg_vals[data$from] + deg_vals[data$to]
      degree_contrib[is.na(degree_contrib)] <- 0

      deg_contribs[["degree_effects"]] <- degree_contrib
      lp_degree <- degree_contrib
    }
  }

  # 3. Compute Baseline/Temporal effects contributions
  baseline_contrib <- rep(0, nrow(data))
  est_time <- model$est_time
  time_changepoints <- model$time_changepoints

  if (!is.null(est_time) && length(est_time) > 0 && !is.null(time_changepoints)) {
    gamma <- if (isTRUE(model$full_baseline)) est_time else c(0, est_time)
    bins <- findInterval(data$time, c(-Inf, time_changepoints, Inf))
    contrib_block <- gamma[bins]
    contrib_block[is.na(contrib_block)] <- 0
    baseline_contrib <- baseline_contrib + contrib_block
  }
  if (length(is_time) > 0) {
    existing_time <- intersect(nm[is_time], names(data))
    if (length(existing_time) > 0) {
      X_time <- as.matrix(data[, existing_time, with = FALSE])
      beta_time <- model$coefficients[existing_time]
      contrib_param <- as.vector(X_time %*% beta_time)
      contrib_param[is.na(contrib_param)] <- 0
      baseline_contrib <- baseline_contrib + contrib_param
    }
  }

  names_vector <- if (!is.null(data$from) && !is.null(data$to)) {
    paste(data$from, data$to, sep = "->")
  } else {
    NULL
  }

  # 4. Assemble linear predictor
  lp <- lp_core + lp_degree + baseline_contrib

  if (type == "lp") {
    names(lp) <- names_vector
    return(lp)
  } else if (type == "response") {
    res <- exp(lp)
    names(res) <- names_vector
    return(res)
  } else if (type == "terms") {
    # Combine all contributions into a matrix/data frame
    res_list <- list()
    # Add core terms
    for (n in names(core_contribs)) {
      res_list[[n]] <- core_contribs[[n]]
    }
    # Add degree effects
    for (n in names(deg_contribs)) {
      res_list[[n]] <- deg_contribs[[n]]
    }
    # Add time effects if present in model
    if ((!is.null(est_time) && length(est_time) > 0 && !is.null(time_changepoints)) || length(is_time) > 0) {
      res_list[["time_effects"]] <- baseline_contrib
    }

    res_mat <- do.call(cbind, res_list)
    if (is.null(res_mat)) {
      res_mat <- matrix(0, nrow = nrow(data), ncol = 0)
    }
    rownames(res_mat) <- names_vector
    return(res_mat)
  }
}

#' The dem Object
#'
#' An object of class \code{dem} returned by the \code{\link{dem}} function,
#' representing a fitted Durational Event Model.
#'
#' @section Value:
#' A \code{dem} object is a list containing the following components:
#' \itemize{
#'   \item \code{call}: The matched call.
#'   \item \code{event_numbers}: A vector containing the number of events.
#'   \item \code{model_0_1}: The fitted model for transition 0 -> 1 (formation).
#'   \item \code{model_1_0}: The fitted model for transition 1 -> 0 (dissolution).
#'   \item \code{events}: The preprocessed event matrix.
#'   \item \code{formula_0_1}: The formula for transition 0 -> 1.
#'   \item \code{formula_1_0}: The formula for transition 1 -> 0.
#'   \item \code{n_nodes}: The number of nodes.
#'   \item \code{simultaneous_interactions}: Logical indicating whether
#'     simultaneous interactions were allowed.
#'   \item \code{directed}: Logical indicating whether the events are directed.
#'   \item \code{training_start}: The start time of the training period.
#'   \item \code{build_time}: The time at which the estimation dataset started
#'     building.
#'   \item \code{max_time}: The maximum event time.
#'   \item \code{exogenous_end}: The end time of the exogenous period.
#'   \item \code{time_changepoints}: Time points where baseline intensity changes.
#'   \item \code{labels_changepoints}: Labels for the time intervals.
#'   \item \code{subsample}: Subsample proportion used.
#'   \item \code{return_data}: Logical indicating whether preprocessed data
#'     frames were returned.
#'   \item \code{runtime}: The estimation runtime.
#'   \item \code{window_map}: The window map used for calculation.
#'   \item \code{preprocessed}: Preprocessed data structures.
#' }
#'
#' @section Methods (S3):
#' The following S3 methods are implemented for \code{dem} objects:
#' \itemize{
#'   \item \code{print(x, ...)}: Prints a brief summary of the DEM model.
#'     \itemize{
#'       \item \code{x}: A \code{dem} object.
#'       \item \code{...}: Additional arguments passed to printing function.
#'     }
#'   \item \code{summary(object, ...)}: Summarizes model results, including
#'     parameter estimates, standard errors, and fit statistics (AIC/BIC).
#'     Returns an object of class \code{summary.dem}.
#'     \itemize{
#'       \item \code{object}: A \code{dem} object.
#'       \item \code{...}: Additional arguments passed to summary method.
#'     }
#'   \item \code{plot(x, which = 3, separate = FALSE, baseline = FALSE, ...)}:
#'     Generates trace plots for coefficients and log-likelihood histories.
#'     \itemize{
#'       \item \code{x}: A \code{dem} object.
#'       \item \code{which}: Integer indicating transition plots to display.
#'         Options are \code{1} for formation (0 -> 1), \code{2} for
#'         dissolution (1 -> 0), or \code{3} (default) for both.
#'       \item \code{separate}: Logical. If \code{TRUE}, each plot is
#'         generated in a new window or as a separate sequence.
#'         Defaults to \code{FALSE}.
#'       \item \code{baseline}: Logical. If \code{TRUE}, plots the estimated
#'         baseline step function. Defaults to \code{FALSE}.
#'       \item \code{...}: Additional arguments passed to the underlying
#'         trace plotting method. Supported parameters include:
#'         \itemize{
#'           \item \code{coefs}: Logical. If \code{TRUE} (default), plots
#'             iteration traces for core/fixed coefficients.
#'           \item \code{degree}: Logical. If \code{TRUE} (default), plots
#'             iteration traces for degree/actor effects.
#'           \item \code{time}: Logical. If \code{TRUE} (default), plots
#'             iteration traces for temporal baseline effects.
#'           \item \code{llh}: Logical. If \code{TRUE} (default), plots
#'             trace of log-likelihood history.
#'           \item \code{sub_label}: Character. Optional subtitle to display
#'             at the bottom of the figure.
#'         }
#'     }
#'   \item \code{plot_baseline(x, process = c("formation", "dissolution"), ...)}:
#'     Plots the step function of the estimated baseline intensity.
#'     \itemize{
#'       \item \code{x}: A \code{dem} object.
#'       \item \code{process}: Character. Specifies whether to plot the
#'         baseline for \code{"formation"} (0 -> 1) (default) or the
#'         \code{"dissolution"} (1 -> 0) process.
#'       \item \code{...}: Additional graphical parameters passed to
#'         \code{\link[graphics]{plot}}.
#'     }
#'   \item \code{predict(object, time = NULL, type = c("response", "lp", "terms"), process = c("both", "formation", "dissolution"), ...)}:
#'     Predicts the intensity, linear predictor, or term contributions for a fitted DEM model.
#'     \itemize{
#'       \item \code{object}: A \code{dem} object.
#'       \item \code{time}: Numeric vector; optional time point(s) at which to predict. Defaults to NULL.
#'       \item \code{type}: Character; the type of prediction. Defaults to \code{"response"}.
#'       \item \code{process}: Character; the transition process to predict. Defaults to \code{"both"}.
#'       \item \code{...}: Additional arguments.
#'     }
#' }
#'
#' @name dem_object
#' @aliases dem-class
NULL

#' The rem Object
#'
#' An object of class \code{rem} returned by the \code{\link{rem}} function,
#' representing a fitted Relational Event Model.
#'
#' @section Value:
#' A \code{rem} object is a list containing the following components:
#' \itemize{
#'   \item \code{call}: The matched call.
#'   \item \code{event_numbers}: A vector containing the number of events.
#'   \item \code{model}: The fitted underlying model.
#'   \item \code{events}: The preprocessed event matrix.
#'   \item \code{formula}: The formula used.
#'   \item \code{n_nodes}: The number of nodes.
#'   \item \code{directed}: Logical indicating whether the events are directed.
#'   \item \code{build_time}: The time at which the estimation dataset started
#'     building.
#'   \item \code{max_time}: The maximum event time.
#'   \item \code{time_changepoints}: Time points where baseline intensity changes.
#'   \item \code{labels_changepoints}: Labels for the time intervals.
#'   \item \code{training_start}: The start time of the training period.
#'   \item \code{exogenous_end}: The end time of the exogenous period.
#'   \item \code{subsample}: Subsample proportion used.
#'   \item \code{return_data}: Logical indicating whether preprocessed data
#'     frames were returned.
#'   \item \code{runtime}: The estimation runtime.
#'   \item \code{window_map}: The window map used for calculation.
#'   \item \code{preprocessed}: Preprocessed data structures.
#' }
#'
#' @section Methods (S3):
#' The following S3 methods are implemented for \code{rem} objects:
#' \itemize{
#'   \item \code{print(x, ...)}: Prints a brief summary of the fitted REM model.
#'     \itemize{
#'       \item \code{x}: A \code{rem} object.
#'       \item \code{...}: Additional arguments passed to printing function.
#'     }
#'   \item \code{summary(object, ...)}: Summarizes model results, including
#'     parameter estimates, standard errors, and fit statistics. Returns
#'     an object of class \code{summary.rem}.
#'     \itemize{
#'       \item \code{object}: A \code{rem} object.
#'       \item \code{...}: Additional arguments passed to summary method.
#'     }
#'   \item \code{plot(x, baseline = FALSE, ...)}: Generates trace plots for
#'     model coefficients and log-likelihood histories.
#'     \itemize{
#'       \item \code{x}: A \code{rem} object.
#'       \item \code{baseline}: Logical. If \code{TRUE}, plots the estimated
#'         baseline step function. Defaults to \code{FALSE}.
#'       \item \code{...}: Additional arguments passed to underlying trace
#'         plotting method. Supported parameters include:
#'         \itemize{
#'           \item \code{coefs}: Logical. If \code{TRUE} (default), plots
#'             iteration traces for core/fixed coefficients.
#'           \item \code{degree}: Logical. If \code{TRUE} (default), plots
#'             iteration traces for degree/actor effects.
#'           \item \code{time}: Logical. If \code{TRUE} (default), plots
#'             iteration traces for temporal baseline effects.
#'           \item \code{llh}: Logical. If \code{TRUE} (default), plots
#'             trace of log-likelihood history.
#'           \item \code{separate}: Logical. If \code{TRUE}, each plot is
#'             generated in a new window or as a separate sequence.
#'             Defaults to \code{FALSE}.
#'           \item \code{sub_label}: Character. Optional subtitle to display
#'             at the bottom of the figure.
#'         }
#'     }
#'   \item \code{logLik(object, ...)}: Extracts the log-likelihood of the
#'     fitted model.
#'     \itemize{
#'       \item \code{object}: A \code{rem} object.
#'       \item \code{...}: Additional arguments.
#'     }
#'   \item \code{plot_baseline(x, ...)}: Plots the step function of the
#'     estimated baseline intensity.
#'     \itemize{
#'       \item \code{x}: A \code{rem} object.
#'       \item \code{...}: Additional graphical parameters passed to
#'         \code{\link[graphics]{plot}}.
#'     }
#'   \item \code{predict(object, time = NULL, type = c("response", "lp", "terms"), ...)}:
#'     Predicts the intensity, linear predictor, or term contributions for a fitted REM model.
#'     \itemize{
#'       \item \code{object}: A \code{rem} object.
#'       \item \code{time}: Numeric vector; optional time point(s) at which to predict. Defaults to NULL.
#'       \item \code{type}: Character; the type of prediction. Defaults to \code{"response"}.
#'       \item \code{...}: Additional arguments.
#'     }
#' }
#'
#' @name rem_object
#' @aliases rem-class
NULL
