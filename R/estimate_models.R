# Helper for safe maximum calculation
safe_max <- function(x, default = 0) {
  x <- x[is.finite(x)]
  if (length(x) > 0) max(x) else default
}

identify_identifiable <- function(y, X) {
  # 1. Validation
  if (length(y) != nrow(X)) {
    stop("Dimensions mismatch: length(y) must equal nrow(X).")
  }
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  # 2. Define subset for events
  idx_pos <- y > 0

  # 3. Compute indicators for non-identifiability

  # Case: Unilateral Separation (MLE -> -Inf)
  # TRUE if the covariate is always zero when y > 0
  is_separated <- colSums(X[idx_pos, , drop = FALSE] != 0) == 0

  # Case: Global Invariance (Collinearity with Intercept)
  # TRUE if the covariate has no variance across the entire sample
  # Using a single pass for min/max is faster than unique() for large N
  glob_min <- apply(X, 2, min)
  glob_max <- apply(X, 2, max)
  is_constant <- glob_min == glob_max

  # 4. Identify indices that pass both checks
  # Covariates must NOT be separated AND NOT be constant
  identifiable_indices <- which(!(is_separated | is_constant))

  # Ensure the result is a plain integer vector
  return(as.vector(identifiable_indices))
}

#' MM Algorithm for Durational Event Models with Time-Varying Effects
#'
#' Implementation of the scalable block-coordinate ascent algorithm for DEMs
#' with time-varying baseline intensities. This function performs iterative
#' Minorization-Maximization (MM) updates for degree and temporal effects
#' while using Newton-Raphson for core covariates.
#'
#' @details
#' The algorithm decomposes the log-likelihood and updates blocks of parameters
#' sequentially. Specifically:
#' \enumerate{
#'   \item Core effects (\eqn{\beta}) are updated via Newton-Raphson.
#'   \item Degree effects (\eqn{\alpha}) are updated using an MM step that
#'         avoids explicit Hessian inversion for high-dimensional actor sets.
#'   \item Temporal effects (\eqn{\gamma}) are updated via a similar MM step
#'         across defined time changepoints.
#' }
#'
#' @param data Preprocessed data.table.
#' @param indicators Numeric vector of covariate indices.
#' @param it_max Maximum number of iterations.
#' @param n_nodes Number of nodes.
#' @param tol Convergence tolerance.
#' @param accelerated Logical; use SQUAREM acceleration for degree effects.
#' @param time_changepoints Numeric vector of time changepoints.
#' @param labels_changepoints Character vector of labels for time slices.
#' @param subsample Subsampling rate for GLM backup estimation.
#' @param verbose Logical; print progress.
#' @param est_degree Initial degree coefficients.
#' @param est_core Initial core coefficients.
#' @param est_time Initial time effects.
#' @param estimate_degree Logical; estimate degree effects.
#' @param directed Logical; whether the network is directed.
#' @param return_data Logical; whether to return the preprocessed data in the result.
#' @param save_hist Logical; whether to save the iteration history of coefficients.
#' @param use_glm Logical; whether to use GLM-based core updates as fallback or control.
#' @param inf_unidentifiable Logical; if TRUE, unidentifiable parameters are set to -Inf.
#' @keywords internal
estimate_mmt <- function(data,
                         indicators,
                         it_max,
                         n_nodes,
                         tol = 1e-10,
                         accelerated = TRUE,
                         time_changepoints = NULL,
                         labels_changepoints = NULL,
                         subsample = 0.2,
                         verbose = FALSE,
                         est_degree = NULL,
                         est_core = NULL,
                         est_time = NULL,
                         estimate_degree = TRUE,
                         directed = FALSE,
                         return_data = TRUE,
                         save_hist = TRUE,
                         use_glm = FALSE,
                         inf_unidentifiable = TRUE) {
  w <- if ("weight" %in% names(data)) as.numeric(data$weight) else rep(1, nrow(data))
  estimate_time <- !is.null(time_changepoints)
  if (is.null(indicators)) indicators <- character(0)
  indicators <- indicators[!is.na(indicators)]
  covarites <- as.matrix(data[, indicators, with = FALSE])
  if (length(indicators) == 0) covarites <- matrix(numeric(0), nrow = nrow(data), ncol = 0)
  intercept_model <- (ncol(covarites) == 1 && any(tolower(colnames(covarites)) == "intercept"))

  if (estimate_time) {
    has_intercept <- any(tolower(colnames(covarites)) == "intercept")
    full_baseline <- !(estimate_degree || has_intercept)

    if (is.null(labels_changepoints)) {
      labels_changepoints <- as.character(time_changepoints)
    }
    # Prepare the outcome
    time_slices <- cut(
      data$time_new,
      breaks = c(-Inf, time_changepoints, Inf),
      labels = c("Beg", labels_changepoints)
    )
    time_slices_from <- cut(
      data$time,
      breaks = c(-Inf, time_changepoints, Inf),
      labels = c("Beg", labels_changepoints)
    )
    time_names <- levels(time_slices)
    if (!full_baseline) time_names <- time_names[-1]
    time_slices_from <- as.numeric(time_slices_from)
    time_slices <- as.numeric(time_slices)

    obs_w <- data$event * (if ("weight" %in% names(data)) data$weight else 1)
    weight_w <- (if ("weight" %in% names(data)) data$weight else 1)
    update_data_time <- get_data_time_alt(
      slice_start_v_r = time_slices_from,
      slice_end_v_r = time_slices,
      time_start_v_r = data$time,
      time_end_v_r = data$time_new,
      weight_v_r = rep(0, nrow(data)), # Initial pred weight is 0
      observation_v_r = obs_w,
      n_slices = length(time_changepoints) + 1,
      changepoints_r = c(time_changepoints),
      full_baseline = full_baseline
    )

    # Check if we have any temporal parameters to validate
    while (nrow(update_data_time) > 0 && ncol(update_data_time) >= 2) {
      is_identifiable_time <- update_data_time[, 2] != 0
      if (sum(!is_identifiable_time) == 0) break

      to_remove_idx <- which(!is_identifiable_time)
      cp_to_remove <- unique(pmin(pmax(to_remove_idx, 1), length(time_changepoints)))

      if (length(cp_to_remove) > 0) {
        labels_changepoints <- labels_changepoints[-cp_to_remove]
        time_changepoints <- time_changepoints[-cp_to_remove]
        if (!is.null(est_time) && length(est_time) > 0) {
          if (full_baseline) {
            est_time <- est_time[-to_remove_idx]
          } else {
            est_time <- est_time[-pmax(1, to_remove_idx - 1)]
          }
        }
      }

      breaks <- c(-Inf, time_changepoints, Inf)
      time_slices <- cut(data$time_new, breaks = breaks, labels = c("Beg", labels_changepoints))
      time_slices_from <- cut(data$time, breaks = breaks, labels = c("Beg", labels_changepoints))
      time_names <- levels(time_slices)
      if (!full_baseline) time_names <- time_names[-1]
      time_slices_from <- as.numeric(time_slices_from)
      time_slices <- as.numeric(time_slices)

      obs_w <- data$event * (if ("weight" %in% names(data)) data$weight else 1)
      update_data_time <- get_data_time_alt(
        slice_start_v_r = time_slices_from,
        slice_end_v_r = time_slices,
        time_start_v_r = data$time,
        time_end_v_r = data$time_new,
        weight_v_r = rep(0, nrow(data)),
        observation_v_r = obs_w,
        n_slices = length(time_changepoints) + 1,
        changepoints_r = c(time_changepoints),
        full_baseline = full_baseline
      )
    }
    data$time_slices_from <- time_slices_from
    data$time_slices <- time_slices
    if (is.null(est_time) || length(est_time) == 0) {
      est_time <- numeric(if (full_baseline) length(time_changepoints) + 1 else length(time_changepoints))
    }

    if (length(est_time) > 0) {
      offset_time <- log(get_time_offset(
        from_slice_r = time_slices_from,
        to_slice_r = time_slices,
        from_time_r = data$time,
        to_time_r = data$time_new,
        est_time_r = if (full_baseline) exp(est_time) else exp(c(0, est_time)),
        changepoints_r = c(time_changepoints, safe_max(data$time_new))
      ))
    } else {
      offset_time <- data$offset
    }
  } else {
    full_baseline <- FALSE
    offset_time <- data$offset
    time_names <- NULL
    est_time_hist <- NULL
  }

  # 1. Initialize parameters and names
  if (is.null(est_core) || length(est_core) == 0) {
    est_core <- rep(0, length(indicators))
  }
  est_core_names <- if (is.character(indicators)) indicators else names(data)[indicators]
  if (!is.null(est_core_names) && length(est_core) == length(est_core_names)) {
    names(est_core) <- est_core_names
  }

  n_core <- length(est_core)

  if (is.null(est_degree) || length(est_degree) == 0) {
    est_degree <- rep(if (estimate_degree) -5 else 0, if (directed) 2 * n_nodes else n_nodes)
  }

  if (intercept_model) {
    est_core <- 0
    names(est_core) <- "Intercept"
  }

  if (directed) {
    if (length(est_degree) == n_nodes) {
      est_mu <- est_degree
      est_nu <- est_degree
    } else {
      est_mu <- est_degree[1:n_nodes]
      est_nu <- est_degree[(n_nodes + 1):(2 * n_nodes)]
    }
  }

  est_degree_names <- if (directed) {
    c(paste0("sender_", seq_len(n_nodes)), paste0("receiver_", seq_len(n_nodes)))
  } else {
    paste0("effect_", seq_len(n_nodes))
  }
  names(est_degree) <- est_degree_names

  if (estimate_time && (is.null(est_time) || length(est_time) == 0)) {
    est_time <- rep(0, if (full_baseline) length(time_changepoints) + 1 else length(time_changepoints))
  }
  if (estimate_time) {
    if (length(est_time) != length(time_names)) {
      # warning(paste("Length mismatch for est_time:", length(est_time), "vs time_names:", length(time_names)))
      if (length(est_time) < length(time_names)) {
        est_time <- c(est_time, rep(0, length(time_names) - length(est_time)))
      } else {
        est_time <- est_time[1:length(time_names)]
      }
    }
    names(est_time) <- paste0("time_", time_names)
  }

  # 2. Initialize history matrices
  if (save_hist) {
    est_degree_hist <- matrix(0, nrow = it_max, ncol = length(est_degree))
    colnames(est_degree_hist) <- est_degree_names

    coefficients_core_hist <- matrix(0, nrow = it_max, ncol = length(est_core))
    colnames(coefficients_core_hist) <- est_core_names

    if (estimate_time) {
      est_time_hist <- matrix(0, nrow = it_max, ncol = length(est_time))
      colnames(est_time_hist) <- names(est_time)
    }
  } else {
    est_degree_hist <- coefficients_core_hist <- est_time_hist <- NULL
  }

  est_mu_old <- if (directed) est_mu else NULL
  est_nu_old <- if (directed) est_nu else NULL
  est_degree_old <- if (!directed) est_degree else NULL

  if (intercept_model) {
    offset_core <- rep(as.vector(est_core), nrow(data))
  } else {
    offset_core <- as.vector(covarites %*% est_core)
  }
  if (length(offset_core) == 0) offset_core <- rep(0, nrow(data))

  offset_degree <- if (directed) as.vector(est_mu[data$from] + est_nu[data$to]) else as.vector(est_degree[data$from] + est_degree[data$to])
  if (length(offset_degree) == 0) offset_degree <- rep(0, nrow(data))

  if (estimate_time) {
    if(all(est_time == 0)){
      offset_time <- rep(0, nrow(data))
    } else {
      offset_time <- log(get_time_offset(
        from_slice_r = data$time_slices_from,
        to_slice_r = data$time_slices,
        from_time_r = data$time,
        to_time_r = data$time_new,
        est_time_r = if (full_baseline) exp(est_time) else exp(c(0, est_time)),
        changepoints_r = c(time_changepoints, safe_max(data$time_new))
      ))
    }
  } else {
    offset_time <- if (!is.null(data$offset)) as.vector(data$offset) else rep(0, nrow(data))
  }


  # if (length(offset_time) == 0) offset_time <- rep(0, nrow(data))

  llh_hist <- numeric(it_max)

  prediction <- exp(as.vector(offset_time) + as.vector(offset_degree) + as.vector(offset_core))
  prediction[!is.finite(prediction)] <- 0

  if (ncol(covarites) > 0) {
    identifiable <- logical(ncol(covarites))
    identifiable[identify_identifiable(y = data$event, X = covarites)] <- TRUE
  } else {
    identifiable <- logical(0)
  }

  for (it in 1:it_max) {
    llh_before <- eval_llh_pois(outcome = data$event, mean = prediction, weights = w)
    if (verbose) {
      cat("\rIteration:", it, "/", it_max, "- LLH:", round(llh_before, 4))
      utils::flush.console()
    }
    if (n_core > 0) {
      offset_fixed <- offset_degree + offset_time

      if (use_glm) {
        est_core <- update_core_glm(
          data = data,
          covarites = covarites,
          est_core = est_core,
          identifiable = identifiable,
          offset_fixed = offset_fixed,
          subsample = subsample
        )
      } else {
        est_core_new <- update_core_cpp(
          X = covarites,
          y = data$event,
          prediction = prediction,
          est_core = est_core,
          identifiable = which(identifiable) - 1,
          offset_fixed = offset_fixed,
          weights = w
        )

        # Check if NR improved likelihood (Reference Fallback)
        prediction_new <- exp(offset_degree + as.vector(covarites %*% est_core_new) + offset_time)
        # prediction_new[!is.finite(prediction_new)] <- safe_max(prediction_new[is.finite(prediction_new)])
        llh_after <- eval_llh_pois(outcome = data$event, mean = prediction_new, weights = w)

        if (!is.finite(llh_after) || llh_after < llh_before - 1e-9) {
          est_core <- update_core_glm(
            data = data,
            covarites = covarites,
            est_core = est_core,
            identifiable = identifiable,
            offset_fixed = offset_fixed,
            subsample = subsample
          )
        } else {
          est_core <- est_core_new
        }
      }
      prediction <- exp(as.vector(offset_degree) + as.vector(covarites %*% est_core) + as.vector(offset_time))
      prediction[!is.finite(prediction)] <- safe_max(prediction[is.finite(prediction)])
      }

    if (save_hist) {
      coefficients_core_hist[it, ] <- est_core
    }
    if (intercept_model) {
      offset_core <- rep(as.vector(est_core), nrow(data))
    } else {
      offset_core <- as.vector(covarites %*% est_core)
    }
    if (length(offset_core) == 0) offset_core <- rep(0, nrow(data))

    if (estimate_degree) {
      # Robustness check for infinite core effects
      offset_fixed <- offset_core + offset_time
      offset_fixed[is.infinite(offset_fixed)] <- safe_max(offset_fixed[is.finite(offset_fixed)])

      if (accelerated) {
        if (directed) {
          # Step 1: x1 = G(x0)
          est_pos <- update_degree_fast(
            from_v = data$from, to_v = data$to,
            event_v = data$event, prediction_v = prediction,
            weights = if ("weight" %in% names(data)) data$weight else numeric(0),
            est_degree = c(est_mu, est_nu), n_nodes = n_nodes,
            directed = TRUE, update_sender = TRUE
          )$est_degree
          
          # We need to update prediction after each coordinate step to be mathematically identical to before
          offset_degree_tmp <- est_pos[data$from] + est_pos[n_nodes + data$to]
          pred_tmp <- exp(offset_degree_tmp + offset_core + offset_time)
          
          est_pos <- update_degree_fast(
            from_v = data$from, to_v = data$to,
            event_v = data$event, prediction_v = pred_tmp,
            weights = if ("weight" %in% names(data)) data$weight else numeric(0),
            est_degree = est_pos, n_nodes = n_nodes,
            directed = TRUE, update_sender = FALSE
          )$est_degree
          
          mu_pos <- est_pos[seq_len(n_nodes)]
          nu_pos <- est_pos[(n_nodes + 1):(2 * n_nodes)]

          # Step 2: x2 = G(x1)
          offset_degree_tmp <- mu_pos[data$from] + nu_pos[data$to]
          pred_tmp <- exp(offset_degree_tmp + offset_core + offset_time)
          
          est_pos_pos <- update_degree_fast(
            from_v = data$from, to_v = data$to,
            event_v = data$event, prediction_v = pred_tmp,
            weights = if ("weight" %in% names(data)) data$weight else numeric(0),
            est_degree = est_pos, n_nodes = n_nodes,
            directed = TRUE, update_sender = TRUE
          )$est_degree
          
          offset_degree_tmp <- est_pos_pos[seq_len(n_nodes)][data$from] + est_pos_pos[(n_nodes + 1):(2*n_nodes)][data$to]
          pred_tmp <- exp(offset_degree_tmp + offset_core + offset_time)
          
          est_pos_pos <- update_degree_fast(
            from_v = data$from, to_v = data$to,
            event_v = data$event, prediction_v = pred_tmp,
            weights = if ("weight" %in% names(data)) data$weight else numeric(0),
            est_degree = est_pos_pos, n_nodes = n_nodes,
            directed = TRUE, update_sender = FALSE
          )$est_degree
          
          mu_pos_pos <- est_pos_pos[seq_len(n_nodes)]
          nu_pos_pos <- est_pos_pos[(n_nodes + 1):(2 * n_nodes)]

          # Acceleration Step: r = x1 - x0, v = (x2 - x1) - r
          r <- est_pos - c(est_mu, est_nu)
          v <- (est_pos_pos - est_pos) - r

          # Scalar alpha calculation with safety epsilon
          alpha <- -sqrt(sum(r^2)) / (sqrt(sum(v^2)) + 1e-10)

          # x_acc = x0 - 2*alpha*r + alpha^2*v
          est_acc_point <- c(est_mu, est_nu) - 2 * alpha * r + alpha^2 * v

          # Final Step: x_final = G(x_acc)
          offset_degree_tmp <- est_acc_point[seq_len(n_nodes)][data$from] + est_acc_point[(n_nodes + 1):(2*n_nodes)][data$to]
          pred_tmp <- exp(offset_degree_tmp + offset_core + offset_time)
          
          est_acc_final <- update_degree_fast(
            from_v = data$from, to_v = data$to,
            event_v = data$event, prediction_v = pred_tmp,
            weights = if ("weight" %in% names(data)) data$weight else numeric(0),
            est_degree = est_acc_point, n_nodes = n_nodes,
            directed = TRUE, update_sender = TRUE
          )$est_degree
          
          offset_degree_tmp <- est_acc_final[seq_len(n_nodes)][data$from] + est_acc_final[(n_nodes + 1):(2*n_nodes)][data$to]
          pred_tmp <- exp(offset_degree_tmp + offset_core + offset_time)
          
          est_acc_final <- update_degree_fast(
            from_v = data$from, to_v = data$to,
            event_v = data$event, prediction_v = pred_tmp,
            weights = if ("weight" %in% names(data)) data$weight else numeric(0),
            est_degree = est_acc_final, n_nodes = n_nodes,
            directed = TRUE, update_sender = FALSE
          )$est_degree

          # Likelihood check: Compare against the second MM step (x2)
          llh_acc <- eval_llh_pois_log(data$event, est_acc_final[data$from] + est_acc_final[data$to + n_nodes] + offset_core + offset_time, weights = w)
          llh_mm <- eval_llh_pois_log(data$event, mu_pos_pos[data$from] + nu_pos_pos[data$to] + offset_core + offset_time, weights = w)

          if (llh_acc > llh_mm) {
            est_degree <- est_acc_final
          } else {
            est_degree <- est_pos_pos
          }

          est_mu <- est_degree[seq_len(n_nodes)]
          est_nu <- est_degree[(n_nodes + 1):(2 * n_nodes)]

          # Identifiability shift
          finite_nu <- which(is.finite(est_nu) & est_nu > -25)
          if (length(finite_nu) > 0) {
            shift <- est_nu[finite_nu[1]]
            est_nu <- est_nu - shift
            est_mu <- est_mu + shift
          }
          est_degree <- c(est_mu, est_nu)
          offset_degree <- est_mu[data$from] + est_nu[data$to]
        } else {
          # Accelerated update (SQUAREM)
          # Step 1: x1 = G(x0)
          est_pos <- update_degree_fast(
            from_v = data$from, to_v = data$to,
            event_v = data$event, prediction_v = prediction,
            weights = if ("weight" %in% names(data)) data$weight else numeric(0),
            est_degree = est_degree, n_nodes = n_nodes,
            directed = FALSE
          )$est_degree

          # Step 2: x2 = G(x1)
          offset_degree_tmp <- est_pos[data$from] + est_pos[data$to]
          pred_tmp <- exp(offset_degree_tmp + offset_core + offset_time)
          
          est_pos_pos <- update_degree_fast(
            from_v = data$from, to_v = data$to,
            event_v = data$event, prediction_v = pred_tmp,
            weights = if ("weight" %in% names(data)) data$weight else numeric(0),
            est_degree = est_pos, n_nodes = n_nodes,
            directed = FALSE
          )$est_degree

          # Acceleration Step
          r <- est_pos - est_degree
          v <- (est_pos_pos - est_pos) - r
          alpha <- -sqrt(sum(r^2)) / (sqrt(sum(v^2)) + 1e-10)

          # x_acc = x0 - 2*alpha*r + alpha^2*v
          est_acc_point <- est_degree - 2 * alpha * r + alpha^2 * v

          # Final Step: x_final = G(x_acc)
          offset_degree_tmp <- est_acc_point[data$from] + est_acc_point[data$to]
          pred_tmp <- exp(offset_degree_tmp + offset_core + offset_time)
          
          est_acc_final <- update_degree_fast(
            from_v = data$from, to_v = data$to,
            event_v = data$event, prediction_v = pred_tmp,
            weights = if ("weight" %in% names(data)) data$weight else numeric(0),
            est_degree = est_acc_point, n_nodes = n_nodes,
            directed = FALSE
          )$est_degree

          # Likelihood check: Compare against second MM step (x2)
          offset_degree_acc <- est_acc_final[data$from] + est_acc_final[data$to]
          llh_acc <- eval_llh_pois_log(outcome = data$event, log_mean = offset_degree_acc + offset_core + offset_time, weights = w)

          offset_degree_mm <- est_pos_pos[data$from] + est_pos_pos[data$to]
          llh_mm <- eval_llh_pois_log(outcome = data$event, log_mean = offset_degree_mm + offset_core + offset_time, weights = w)

          if (llh_acc > llh_mm) {
            est_degree <- est_acc_final
          } else {
            est_degree <- est_pos_pos
          }
          offset_degree <- est_degree[data$from] + est_degree[data$to]
        }
      } else {
        # Standard MM update
        if (directed) {
          # Update sender then receiver
          est_degree <- update_degree_fast(
            from_v = data$from, to_v = data$to,
            event_v = data$event, prediction_v = prediction,
            weights = if ("weight" %in% names(data)) data$weight else numeric(0),
            est_degree = c(est_mu, est_nu), n_nodes = n_nodes,
            directed = TRUE, update_sender = TRUE
          )$est_degree
          
          offset_degree_tmp <- est_degree[seq_len(n_nodes)][data$from] + est_degree[(n_nodes + 1):(2*n_nodes)][data$to]
          pred_tmp <- exp(offset_degree_tmp + offset_core + offset_time)
          
          est_degree <- update_degree_fast(
            from_v = data$from, to_v = data$to,
            event_v = data$event, prediction_v = pred_tmp,
            weights = if ("weight" %in% names(data)) data$weight else numeric(0),
            est_degree = est_degree, n_nodes = n_nodes,
            directed = TRUE, update_sender = FALSE
          )$est_degree
          
          est_mu <- est_degree[seq_len(n_nodes)]
          est_nu <- est_degree[(n_nodes + 1):(2 * n_nodes)]
          
          # Identifiability shift
          finite_nu <- which(is.finite(est_nu) & est_nu > -25)
          if (length(finite_nu) > 0) {
            shift <- est_nu[finite_nu[1]]
            est_nu <- est_nu - shift
            est_mu <- est_mu + shift
          }
          est_degree <- c(est_mu, est_nu)
          offset_degree <- est_mu[data$from] + est_nu[data$to]
        } else {
          est_degree <- update_degree_fast(
            from_v = data$from, to_v = data$to,
            event_v = data$event, prediction_v = prediction,
            weights = if ("weight" %in% names(data)) data$weight else numeric(0),
            est_degree = est_degree, n_nodes = n_nodes,
            directed = FALSE
          )$est_degree
          
          est_degree[est_degree < -50] <- -50
          offset_degree <- est_degree[data$from] + est_degree[data$to]
        }
      }
    }
    if (length(offset_degree) == 0) offset_degree <- rep(0, nrow(data))
    if (save_hist) {
      est_degree_hist[it, ] <- est_degree
    }

    if (estimate_time) {
      # Robustness for time update
      pred_pop_core <- exp(offset_core + offset_degree)
      pred_pop_core[is.infinite(pred_pop_core)] <- safe_max(pred_pop_core[is.finite(pred_pop_core)])

      obs_w <- data$event * (if ("weight" %in% names(data)) data$weight else 1)
      pred_w <- pred_pop_core * (if ("weight" %in% names(data)) data$weight else 1)
      update_data_time <- get_data_time_alt(
        slice_start_v_r = data$time_slices_from,
        slice_end_v_r = data$time_slices,
        time_start_v_r = data$time,
        time_end_v_r = data$time_new,
        weight_v_r = pred_w,
        observation_v_r = obs_w,
        n_slices = length(time_changepoints) + 1,
        changepoints_r = c(time_changepoints),
        full_baseline = full_baseline
      )

      if (nrow(update_data_time) > 0 && ncol(update_data_time) >= 3) {
        # Avoid log(0) or 0/0
        obs <- update_data_time[, 2]
        pred <- update_data_time[, 3]
        pred[pred < 1e-10] <- 1e-10
        est_time <- log(obs / pred)
        # Use a lower bound during iterations
        est_time[est_time < -50] <- -50
      } else {
        est_time <- rep(-10, if (full_baseline) length(time_changepoints) + 1 else length(time_changepoints))
      }

      offset_time <- log(get_time_offset(
        from_slice_r = data$time_slices_from,
        to_slice_r = data$time_slices,
        from_time_r = data$time,
        to_time_r = data$time_new,
        est_time_r = if (full_baseline) exp(est_time) else exp(c(0, est_time)),
        changepoints_r = c(time_changepoints, safe_max(data$time_new))
      ))
      if (save_hist) {
        est_time_hist[it, ] <- est_time
      }
    }
    prediction <- exp(offset_degree + offset_core + offset_time)
    prediction[!is.finite(prediction)] <- safe_max(prediction[is.finite(prediction)])
    
    if (estimate_time) {
      gamma_full <- if (full_baseline) est_time else c(0, est_time)
      event_idx <- which(data$event > 0)
      if (length(event_idx) > 0) {
        slice_idx <- as.integer(data$time_slices[event_idx])
        inst_log_rate <- (if (length(offset_degree) > 1) offset_degree[event_idx] else offset_degree) +
                         (if (length(offset_core) > 1) offset_core[event_idx] else offset_core) +
                         gamma_full[slice_idx]
        llh_iter_end <- sum(w[event_idx] * data$event[event_idx] * inst_log_rate) - sum(w * prediction)
      } else {
        llh_iter_end <- -sum(w * prediction)
      }
    } else {
      llh_iter_end <- eval_llh_pois(outcome = data$event, mean = prediction, weights = w) - sum(w * data$event * data$offset, na.rm = TRUE)
    }
    llh_hist[it] <- llh_iter_end

    if (it > 1) {
      part_1 <- sqrt(sum((est_core - est_core_old)^2) + sum((est_degree - est_degree_old)^2) + (if (estimate_time) sum((est_time - est_time_old)^2) else 0))
      part_2 <- abs(llh_hist[it] - llh_hist[it - 1]) / (abs(llh_hist[it]) + 1e-6)
      val_check <- max(part_1, part_2)
      if (is.na(val_check)) val_check <- Inf
      if (verbose) {
        cat(" - Criterion:", round(val_check, 4))
        utils::flush.console()
      }
      if (val_check < tol) {
        if (verbose) cat("\nConverged after", it, "iterations\n")
        break
      }
    }
    est_core_old <- est_core
    est_degree_old <- est_degree
    if (estimate_time) est_time_old <- est_time
  }
  if (verbose) cat("\n")

  if (estimate_degree) {
    # Final pass: Set unidentifiable degree estimates to -Inf
    obs_counts <- if (directed) {
      c(
        tabulate(data$from[data$event == 1], nbins = n_nodes),
        tabulate(data$to[data$event == 1], nbins = n_nodes)
      )
    } else {
      # For undirected, we sum from and to counts
      counts_from <- tabulate(data$from[data$event == 1], nbins = n_nodes)
      counts_to <- tabulate(data$to[data$event == 1], nbins = n_nodes)
      counts_from + counts_to
    }
    if (inf_unidentifiable) {
      est_degree[obs_counts == 0] <- -Inf
    }

    if (directed) {
      est_mu <- est_degree[1:n_nodes]
      est_nu <- est_degree[(n_nodes + 1):(2 * n_nodes)]
      finite_mu <- which(is.finite(est_mu))
      if (length(finite_mu) > 0) {
        shift <- est_mu[finite_mu[1]]
        est_mu <- est_mu - shift
        est_nu <- est_nu + shift
      }
      est_degree <- c(est_mu, est_nu)
    }
  }

  if (inf_unidentifiable) {
    est_core[!identifiable] <- -Inf
  }
  if (estimate_time) {
    # Final pass for time effects
    pred_degree_core_final <- exp(offset_core + (if (estimate_degree) (if (directed) est_mu[data$from] + est_nu[data$to] else est_degree[data$from] + est_degree[data$to]) else 0))
    obs_w <- data$event * (if ("weight" %in% names(data)) data$weight else 1)
    update_data_time_final <- get_data_time_alt(
      slice_start_v_r = data$time_slices_from, slice_end_v_r = data$time_slices,
      time_start_v_r = data$time, time_end_v_r = data$time_new,
      weight_v_r = pred_degree_core_final, observation_v_r = obs_w,
      n_slices = length(time_changepoints) + 1, changepoints_r = c(time_changepoints),
      full_baseline = full_baseline
    )
    obs_t <- update_data_time_final[, 2]
    if (inf_unidentifiable) {
      # Set unidentifiable time slices to -Inf
      est_time[obs_t == 0] <- -Inf
    }
  }

  # Update prediction one last time before Hessian
  offset_degree <- if (directed) est_degree[data$from] + est_degree[n_nodes + data$to] else est_degree[data$from] + est_degree[data$to]
  if (length(offset_degree) == 0) offset_degree <- rep(0, nrow(data))
  offset_time <- if (estimate_time) {
    log(get_time_offset(
      from_slice_r = data$time_slices_from, to_slice_r = data$time_slices,
      from_time_r = data$time, to_time_r = data$time_new,
      est_time_r = if (full_baseline) exp(est_time) else exp(c(0, est_time)),
      changepoints_r = c(time_changepoints, safe_max(data$time_new))
    ))
  } else {
    if (!is.null(data$offset)) as.vector(data$offset) else rep(0, nrow(data))
  }
  prediction <- exp(offset_degree + offset_core + offset_time)
  prediction[!is.finite(prediction)] <- 0

  time_slices_arg <- if (estimate_time) data$time_slices else rep(1, nrow(data))
  time_slices_arg[is.na(time_slices_arg)] <- 1

  information <- get_A_B_C_D_E_F_exact(
    from_v = data$from,
    to_v = data$to,
    weight_v = prediction,
    covarites = covarites,
    time_slices = time_slices_arg,
    n_slices = if (estimate_time) length(time_changepoints) + 1 else 1,
    n_nodes = n_nodes,
    directed = directed,
    full_baseline = full_baseline
  )

  # Fisher Information and Covariance calculation
  identifiable_pop <- diag(information$B_mat) != 0
  inf_D <- information$D_mat[identifiable_pop, , drop = FALSE]
  inf_F <- information$F_mat[identifiable_pop, , drop = FALSE]
  inf_B <- information$B_mat[identifiable_pop, identifiable_pop, drop = FALSE]

  if (estimate_time) {
    tmp <- sweep(inf_F, 2, 1 / information$C_mat, "*")
    Q_2 <- inf_B - tmp %*% t(inf_F)
    Q_2_inv <- tryCatch(
      {
        if (nrow(Q_2) > 0) solve(Q_2) else matrix(numeric(0), 0, 0)
      },
      error = function(e) {
        if (nrow(Q_2) > 0) MASS::ginv(Q_2) else matrix(numeric(0), 0, 0)
      }
    )
    Y_inv_11 <- t(tmp) %*% Q_2_inv %*% tmp
    diag(Y_inv_11) <- 1 / information$C_mat + diag(Y_inv_11)
    Y_inv_12 <- -t(tmp) %*% Q_2_inv
    fisher_info <- information$A_mat - (information$E_mat %*% Y_inv_11 + t(inf_D) %*% t(Y_inv_12)) %*% t(information$E_mat) -
      (information$E_mat %*% Y_inv_12 + t(inf_D) %*% Q_2_inv) %*% inf_D
  } else {
    A_inv <- tryCatch(
      {
        if (nrow(inf_B) > 0) solve(inf_B) else matrix(numeric(0), 0, 0)
      },
      error = function(e) {
        if (nrow(inf_B) > 0) MASS::ginv(inf_B) else matrix(numeric(0), 0, 0)
      }
    )
    fisher_info <- information$A_mat - t(inf_D) %*% A_inv %*% inf_D
  }

  diag(fisher_info)[diag(fisher_info) == 0] <- 1e-09
  fisher_info_identizable <- fisher_info[identifiable, identifiable, drop = FALSE]
  covariance_core_identifiable <- if (nrow(fisher_info_identizable) > 0) {
    tryCatch(
      {
        solve(fisher_info_identizable)
      },
      error = function(e) MASS::ginv(fisher_info_identizable)
    )
  } else {
    matrix(numeric(0), 0, 0)
  }

  # Final covariance calc
  n_core <- length(est_core)
  covariance_core <- matrix(NA, nrow = n_core, ncol = n_core)
  if (n_core > 0) {
    if (nrow(fisher_info_identizable) > 0) {
      covariance_core[identifiable, identifiable] <- covariance_core_identifiable
    }
  }
  if (!is.null(est_core_names) && length(est_core) == length(est_core_names)) {
    names(est_core) <- est_core_names
  }
  if (!is.null(est_degree_names) && length(est_degree) == length(est_degree_names)) {
    names(est_degree) <- est_degree_names
  }
  if (estimate_time && !is.null(time_names) && length(est_time) == length(time_names)) {
    names(est_time) <- paste0("time_", time_names)
  }

  # Assign column names to history matrices for plotting
  if (save_hist) {
    if (!is.null(est_degree_hist) && !is.null(est_degree_names)) {
      colnames(est_degree_hist) <- est_degree_names
    }
    if (!is.null(coefficients_core_hist) && !is.null(est_core_names)) {
      colnames(coefficients_core_hist) <- est_core_names
    }
    if (estimate_time && !is.null(est_time_hist) && !is.null(time_names)) {
      colnames(est_time_hist) <- paste0("time_", time_names)
    }
  }

  if (data.table::is.data.table(data)) {
    data[, prediction := as.vector(prediction)]
  } else {
    data$prediction <- as.vector(prediction)
  }

  # llh_pois <- eval_llh_pois(outcome = data$event, mean = prediction, weights = w)
  # llh <- llh_pois - sum(w * data$event * data$offset, na.rm = TRUE)
  if (estimate_time) {
    intensity <- exp((if (full_baseline) est_time else c(0, est_time))[data$time_slices] + offset_degree + offset_core)
  } else {
    intensity <- exp(offset_core + offset_degree)
  }


  llh <- calc_llh_scaled(
    pred = data$prediction,
    intensity = intensity,
    delta = data$event,
    pair_id = data$pair_id
  )

  res <- list(
    coefficients = c(est_core, if (estimate_time) est_time else NULL, est_degree),
    est_degree = est_degree,
    est_degree_hist = if (save_hist) est_degree_hist[1:it, , drop = FALSE] else NULL,
    est_core = est_core,
    coefficients_core_hist = if (save_hist) coefficients_core_hist[1:it, , drop = FALSE] else NULL,
    est_time = if (estimate_time) est_time else NULL,
    est_time_hist = if (estimate_time && save_hist) est_time_hist[1:it, , drop = FALSE] else NULL,
    covariance = covariance_core, llh = llh, llh_hist = llh_hist[1:it],
    data = if (return_data) data else NULL,
    n_nodes = n_nodes,
    directed = directed,
    time_changepoints = time_changepoints,
    labels_changepoints = labels_changepoints,
    full_baseline = full_baseline,
    n_obs = nrow(data),
    prediction = as.vector(prediction)
  )
  class(res) <- c("dem.mm", "redeem_result")
  return(res)
}



# estimate_mm <- function(data,
#                         indicators,
#                         it_max,
#                         n_nodes,
#                         tol = 1e-10,
#                         accelerated,
#                         verbose,
#                         subsample,
#                         est_degree,
#                         est_core,
#                         estimate_degree = TRUE,
#                         directed = FALSE,
#                         return_data = TRUE,
#                         save_hist = TRUE) {
#   # Forward call to the unified estimate_mmt engine
#   # When time_changepoints is NULL, estimate_mmt behaves like estimate_mm
#   estimate_mmt(
#     data = data,
#     indicators = indicators,
#     it_max = it_max,
#     n_nodes = n_nodes,
#     tol = tol,
#     accelerated = accelerated,
#     verbose = verbose,
#     subsample = subsample,
#     est_degree = est_degree,
#     est_core = est_core,
#     est_time = NULL, # No temporal baseline coefficients
#     estimate_degree = estimate_degree,
#     directed = directed,
#     time_changepoints = NULL, # No changepoints
#     labels_changepoints = NULL,
#     return_data = return_data,
#     save_hist = save_hist
#   )
# }
#' Update core coefficients using R (for correctness assessment)
#'
#' @inheritParams estimate_mmt
#' @param prediction Current predicted intensities.
#' @param identifiable Logical vector indicating identifiable coefficients.
#' @param offset_fixed Fixed offset (degree + baseline).
#' @keywords internal
update_core_r <- function(data,
                          covarites,
                          prediction,
                          est_core,
                          identifiable,
                          offset_fixed = NULL) {
  if (length(identifiable) == 0 || !any(identifiable)) {
    return(est_core)
  }

  ind_issue <- which(!is.finite(prediction))
  if (length(ind_issue) > 0) {
    prediction[ind_issue] <- safe_max(prediction[-ind_issue])
  }

  # Standard NR update (100% aligned with RILM)
  X_id <- covarites[, identifiable, drop = FALSE]
  tmp_mat <- suppressWarnings(sweep(t(X_id), 2, prediction, "*"))
  H <- tmp_mat %*% X_id

  tmp_mat_inv <- tryCatch(
    {
      solve(H)
    },
    error = function(e) {
      MASS::ginv(H)
    }
  )
  est_core[identifiable] <- est_core[identifiable] +
    as.vector(tmp_mat_inv %*% tmp_mat %*% ((data$event - prediction) / prediction))

  return(est_core)
}

update_core_glm <- function(data,
                            covarites,
                            est_core,
                            identifiable,
                            offset_fixed,
                            subsample = 1.0) {
  if (length(identifiable) == 0 || !any(identifiable)) {
    return(est_core)
  }

  w_val <- if ("weight" %in% names(data)) as.numeric(data$weight) else rep(1, nrow(data))
  n_rows <- nrow(covarites)
  if (subsample < 1.0) {
    subset_idx <- sample.int(n_rows, size = floor(n_rows * subsample))
    y <- (data$event[subset_idx] == 1)
    X_id <- covarites[subset_idx, identifiable, drop = FALSE]
    offset_subset <- offset_fixed[subset_idx]
    w_subset <- w_val[subset_idx]
  } else {
    y <- data$event
    X_id <- covarites[, identifiable, drop = FALSE]
    offset_subset <- offset_fixed
    w_subset <- w_val
  }
  prediction <- exp(offset_fixed + as.vector(covarites %*% est_core))
  llh_before <- eval_llh_pois(outcome = data$event, mean = prediction, weights = w_val)
  suppressWarnings({
    fit <- stats::glm.fit(
      x = X_id,
      y = y,
      offset = offset_subset,
      weights = w_subset,
      family = stats::poisson(),
      control = stats::glm.control(maxit = 1, epsilon = 1e-8)
    )
  })
  new_beta <- est_core
  if (all(is.finite(fit$coefficients))) {
    new_beta[identifiable] <- fit$coefficients
  } else if (any(is.finite(fit$coefficients))) {
    # Match coefficients for identifiable parameters where possible
    idx_finite <- is.finite(fit$coefficients)
    new_beta[identifiable][idx_finite] <- fit$coefficients[idx_finite]
  }
  llh_after <- eval_llh_pois(
    outcome = data$event,
    mean = exp(offset_fixed + as.vector(covarites %*% new_beta)),
    weights = w_val
  )
  if (llh_after < llh_before) {
    suppressWarnings({
      fit <- stats::glm.fit(
        x = X_id,
        y = y,
        offset = offset_subset,
        weights = w_subset,
        family = stats::poisson()
      )
    })

    coef_vals <- fit$coefficients
    coef_vals[is.na(coef_vals)] <- 0
    new_beta <- est_core
    new_beta[identifiable] <- coef_vals
  }
  return(new_beta)
}
