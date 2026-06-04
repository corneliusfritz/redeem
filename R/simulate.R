#                             coef_degree,
#                             time_changepoints,
#                             baseline,
#                             n_nodes,
#                             time,
#                             max_events = 100000,
#                             verbose,
#                             block,
#                             directed = FALSE) {
#   # First, we need to split the between formula into one for each block (TODO)
#   events_between <- rem.simulate(
#     formula = formula_between,
#     coef = coef,
#     coef_degree = coef_degree,
#     time_changepoints = time_changepoints,
#     baseline = baseline,
#     n_nodes = n_nodes,
#     time = time,
#     max_events = max_events,
#     verbose = verbose,
#     block = block,
#     directed = directed
#   )
#   events_within <- list()
#   unique_blocks <- unique(block)
#   # Simulate all the between block events
#   for (i in seq_along(unique_blocks)) {
#     ind <- which(block == unique_blocks[i])
#     # continuous_cov_tmp <- continuous_cov[ind,ind]
#     # categorical_cov_tmp <- categorical_cov[ind,ind]
#     # TODO: Most likely, we need to define a function for within block fomulas and use this here
#     events <- rem.simulate(
#       formula = formula_tmp,
#       coef = coef,
#       coef_popularity = coef_popularity,
#       time_changepoints = time_changepoints,
#       baseline = baseline,
#       n_nodes = length(block_i_nodes),
#       time = time,
#       max_events = max_events,
#       verbose = verbose,
#       directed = directed
#     )
#     events[, 2] <- block_i_nodes[events[, 2]]
#     events[, 3] <- block_i_nodes[events[, 3]]
#     events_within[[i]] <- events
#   }
#   events <- rbind(do.call(rbind, events_within), events_between)
# }
print_simulation_info <- function(preprocessed, coef_0_1_full, coef_1_0_full = NULL,
                                  baseline_0_1, baseline_1_0 = NULL,
                                  all_changepoints, time = 0) {
  cat("\n========================================\n")
  cat("REDEEM Simulation: Terms & Coefficients\n")
  cat("========================================\n")

  # Only show slices up to the simulation time
  if (time > 0) {
    all_changepoints <- all_changepoints[all_changepoints < time]
  }

  slices <- c(0, all_changepoints)
  n_slices <- length(slices)

  for (i in 1:n_slices) {
    t_start <- slices[i]
    t_end <- if (i < n_slices) slices[i + 1] else (if (time > 0) time else Inf)
    cat(sprintf("\n--- Slice %d: [t=%.2f, t=%s] ---\n", i, t_start, if (is.infinite(t_end)) "Inf" else sprintf("%.2f", t_end)))

    # Process 0-1 (Formation in DEM, Event in REM)
    proc_name <- if (is.null(coef_1_0_full)) "Event" else "Formation (0-1)"
    cat(sprintf("  Process: %s\n", proc_name))

    # Use 0-1 specific names if available, otherwise global names (for REM)
    term_names_0_1 <- if (!is.null(preprocessed$preprocess_0_1)) preprocessed$preprocess_0_1$coef_names else preprocessed$coef_names
    for (j in seq_along(term_names_0_1)) {
      # Map process-specific index to global coefficient vector index
      global_idx <- if (!is.null(preprocessed$preprocess_0_1)) {
        match(term_names_0_1[j], preprocessed$coef_names)
      } else {
        j
      }

      if (is.na(global_idx)) next

      val <- coef_0_1_full[global_idx]
      # If it's the intercept, add the baseline shift for this slice
      if (!is.na(term_names_0_1[j]) && tolower(term_names_0_1[j]) == "intercept") {
        val <- val + baseline_0_1[i]
      }
      if (is.na(val) || is.na(term_names_0_1[j])) {
        next
      }
      if (val != 0 || tolower(term_names_0_1[j]) == "intercept") {
        cat(sprintf("    %-20s : %.4f\n", term_names_0_1[j], val))
      }
    }

    # Process 1-0 (Termination in DEM)
    if (!is.null(coef_1_0_full)) {
      cat("\n  Process: Termination (1-0)\n")
      term_names_1_0 <- preprocessed$preprocess_1_0$coef_names
      for (j in seq_along(term_names_1_0)) {
        # Map process-specific index to global coefficient vector index
        global_idx <- match(term_names_1_0[j], preprocessed$coef_names)
        if (is.na(global_idx)) next

        val <- coef_1_0_full[global_idx]
        if (!is.na(term_names_1_0[j]) && tolower(term_names_1_0[j]) == "intercept") {
          val <- val + baseline_1_0[i]
        }
        if (is.na(val) || is.na(term_names_1_0[j])) {
          next
        }
        if (val != 0 || tolower(term_names_1_0[j]) == "intercept") {
          cat(sprintf("    %-20s : %.4f\n", term_names_1_0[j], val))
        }
      }
    }
  }
  cat("========================================\n\n")
}

#' Simulate a Relational Event Model (REM)
#'
#' @details
#' The \code{block} parameter allows the user to specify a partition of the nodes into different groups (blocks).
#' When the vector contains more than one unique block identifier:
#' \itemize{
#'   \item The simulation suppresses all within-block dyad intensities by setting them to 0.
#'   \item Consequently, only events between nodes belonging to different blocks are generated (between-block interactions).
#'   \item If all nodes belong to the same block (e.g., if a single value or \code{NULL} is passed), no block-level constraints are applied, and all dyads are simulated according to the specified model formula.
#' }
#'
#' @param events A matrix representing the initial events with columns \code{time}, \code{from},
#'   \code{to}, and optionally \code{type} (1 for start, 3 for exogenous changes).
#'   Defaults to an empty 4-column matrix.
#' @param formula A one-sided \code{\link[stats]{formula}} specifying the
#'   sufficient statistics to include in the intensity function. The
#'   right-hand side must be composed of terms from
#'   \code{\link{redeem_terms}}. For example:
#'   \code{~ inertia() + reciprocity() + degree}.
#'   An intercept (\code{~ 1}) is the minimal specification.
#' @param coef Numeric vector; coefficients for the model. Defaults to NULL.
#' @param coef_degree Numeric; degree coefficient. Defaults to 0.
#' @param n_events Integer; number of events to simulate. Defaults to 0.
#' @param time Numeric; simulation time. Defaults to 0.
#' @param max_events Integer; maximum number of events. Defaults to 400000.
#' @param n_nodes Integer; the total number of actors in the network.
#' @param verbose Logical; if \code{TRUE}, prints progress information. Defaults to FALSE.
#' @param baseline Numeric vector; baseline intensity values for intervals defined by changepoints. Defaults to NULL.
#' @param seed Integer; random seed. Defaults to 123.
#' @param block An integer vector of length \code{n_nodes} indicating the block/group assignment for each node, or a single value applied to all nodes. Defaults to 1. If multiple blocks are assigned, within-block interactions are suppressed (i.e., their dyadic intensities are set to 0), meaning only events occurring between actors in different blocks are simulated.
#' @param directed Logical; whether the interaction events are directed. Defaults to FALSE.
#' @note Multi-stream event models are currently not supported in simulation.
#' @return A matrix of simulated events.
#' @examples
#' # Simulate events from a REM model structure
#' n <- 10
#' f1 <- ~ 1 + inertia(transformation = "identity")
#'
#' # Simulating events
#' evs <- rem.simulate(
#'   formula = f1,
#'   n_nodes = n,
#'   time = 1.0,
#'   coef = c(1.0, 0.5),
#'   seed = 42,
#'   max_events = 100
#' )
#' head(evs)
#' @export
rem.simulate <- function(events = matrix(0, 0, 4),
                         formula,
                         coef = NULL,
                         coef_degree = 0,
                         n_events = 0,
                         time = 0,
                         max_events = 400000,
                         n_nodes,
                         verbose = FALSE,
                         baseline = NULL,
                         seed = 123, block = 1,
                         directed = FALSE) {
  set.seed(seed)
  if (is.null(coef)) coef <- numeric(0)
  if (is.null(block)) {
    block <- rep(1, n_nodes)
  } else if (length(block) == 1) {
    block <- rep(block, n_nodes)
  }
  if (length(block) != n_nodes) stop("length(block) must be n_nodes")

  preprocessed <- formula_preprocess(
    formula_1_0 = formula,
    formula_0_1 = formula,
    events = events,
    n_nodes = n_nodes,
    directed = directed,
    model_type = "rem",
    simulation = TRUE
  )


  if (any(!sapply(preprocessed$stream_list, is.null))) {
    stop("Multi-stream simulation is not yet supported.")
  }

  if (!isTRUE(preprocessed$preprocess_0_1$includes_degrees)) {
    if (!is.null(coef_degree) && (length(coef_degree) > 1 || any(coef_degree != 0))) {
      stop("coef_degree provided but degrees are not included in the formula.")
    }
  }

  if (length(baseline) > 0 && !isTRUE(preprocessed$preprocess_0_1$has_baseline_term)) {
    stop("baseline estimates were provided but baseline is not included in the formula.")
  }

  coef_augmented <- numeric(length = length(preprocessed$preprocess_0_1$coef_names))
  coef_to_match <- coef
  if (isTRUE(preprocessed$preprocess_0_1$intercept_removed) && length(coef) > 0 && (is.null(names(coef)) || names(coef)[1] == "")) {
    coef_to_match <- coef[-1]
  }
  coef_augmented <- match_coefficients(
    user_coefs = coef_to_match,
    internal_names = preprocessed$preprocess_0_1$coef_names,
    internal_keys = names(preprocessed$preprocess_0_1$coef_names)
  )

  # Map to full vector
  coef_full <- numeric(length(preprocessed$coef_names))
  idx_map <- match(preprocessed$preprocess_0_1$coef_names, preprocessed$coef_names)
  if (any(!is.na(idx_map))) {
    coef_full[idx_map[!is.na(idx_map)]] <- coef_augmented[!is.na(idx_map)]
  }

  # 1. Automatic detection of changepoints from covariates
  cov_changepoints <- NULL
  if (length(preprocessed$data_list) > 0) {
    for (i in seq_along(preprocessed$data_list)) {
      if (is.list(preprocessed$data_list[[i]])) {
        nms <- names(preprocessed$data_list[[i]])
        if (!is.null(nms)) {
          suppressWarnings(ts <- as.numeric(nms))
          ts <- ts[!is.na(ts) & ts > 1e-10]
          cov_changepoints <- c(cov_changepoints, ts)
        }
      }
    }
  }

  all_changepoints <- sort(unique(c(cov_changepoints, preprocessed$baseline_changepoints_0_1)))
  if (is.null(all_changepoints)) all_changepoints <- numeric(0)

  if ((n_events == 0) && (time == 0)) {
    if (length(all_changepoints) > 0) {
      time <- max(all_changepoints)
    } else {
      stop("You need to specify either 'n_events' or 'time'.")
    }
  }

  if ((n_events > 0) && (time > 0)) {
    stop("You need to specify either 'n_events' or 'time', not both.")
  }

  window_info <- if (length(preprocessed$window_map) > 0) {
    as.list(stats::setNames(as.numeric(names(preprocessed$window_map)), preprocessed$window_map))
  } else {
    list()
  }

  if (length(all_changepoints) > 0 || !is.null(baseline)) {
    # Normalize baseline: prepend 0 if not a full baseline (i.e., if intercept/degrees are present)
    has_intercept <- any(tolower(preprocessed$preprocess_0_1$coef_names) == "intercept")
    has_degrees <- isTRUE(preprocessed$preprocess_0_1$includes_degrees) || (!is.null(coef_degree) && any(coef_degree != 0))

    ref_tc <- if (!is.null(preprocessed$baseline_changepoints_0_1)) {
      preprocessed$baseline_changepoints_0_1
    } else {
      numeric(0)
    }

    # If there's no intercept/fixed effects, it's a full baseline.
    full_baseline <- !(has_intercept || has_degrees)

    if (!is.null(baseline)) {
      # If the user provides a full baseline but we have an intercept/degrees,
      # we can potentially shift the first element to the intercept if it's currently 0
      expected_user_len <- if (full_baseline) length(ref_tc) + 1 else length(ref_tc)

      if (!full_baseline && length(baseline) == (length(ref_tc) + 1)) {
        int_idx <- which(tolower(preprocessed$coef_names) == "intercept")
        if (length(int_idx) > 0 && coef_full[int_idx] == 0) {
          coef_full[int_idx] <- baseline[1]
          baseline <- baseline[-1] - baseline[1]
          expected_user_len <- length(ref_tc)
        }
      }

      if (full_baseline && length(baseline) == length(ref_tc)) {
        baseline <- c(0, baseline)
        expected_user_len <- length(ref_tc) + 1
      }

      if (length(baseline) != expected_user_len) {
        stop(sprintf(
          "Expected %d baseline coefficients for %d changepoints (full_baseline=%s), but got %d.",
          expected_user_len, length(ref_tc), full_baseline, length(baseline)
        ))
      }

      baseline_full <- if (full_baseline) baseline else c(0, baseline)

      # Ensure mapping covers the whole time range by adding Inf
      full_tc <- c(-1e-10, ref_tc, Inf)
      slice_indices <- findInterval(c(0, all_changepoints), full_tc)
      baseline <- baseline_full[slice_indices]
    } else {
      baseline <- rep(0, length(all_changepoints) + 1)
    }

    # if (verbose) {
    #   print_simulation_info(
    #     preprocessed = preprocessed,
    #     coef_0_1_full = coef_full,
    #     baseline_0_1 = baseline,
    #     all_changepoints = all_changepoints,
    #     time = time
    #   )
    # }

    res <- rem_simulate_from_empty_timevarying(
      n_nodes = preprocessed$n_nodes,
      coef = coef_full,
      degree_coef = coef_degree,
      baseline = baseline,
      n_events = n_events,
      max_events = max_events,
      time = time,
      time_changepoints = all_changepoints,
      transformations = preprocessed$transformation_list,
      terms = preprocessed$term_names,
      verbose = verbose,
      directed = directed,
      data_list = preprocessed$data_list,
      seed = seed, block = block,
      window_info = window_info
    )
    data_tmp <- res[length(res)]
    res <- res[-length(res)]

    # Calculate intervals
    intervals <- c(0, all_changepoints)
    # Check which time slices have information and exclude the slices with no observations
    res <- do.call(rbind, res)
    return(res)
  } else {
    # if (verbose) {
    #   print_simulation_info(
    #     preprocessed = preprocessed,
    #     coef_0_1_full = coef_full,
    #     baseline_0_1 = 0,
    #     all_changepoints = numeric(0)
    #   )
    # }
    return(
      rem_simulate_from_empty(
        n_nodes = preprocessed$n_nodes,
        coef = coef_full,
        degree_coef = coef_degree,
        n_events = n_events,
        time = time,
        max_events = max_events,
        transformations = preprocessed$transformation_list,
        terms = preprocessed$term_names,
        verbose = verbose,
        directed = directed,
        data_list = preprocessed$data_list,
        seed = seed, block = block,
        window_info = window_info
      )
    )
  }
}

#' Simulate events based on specified formulas and coefficients
#'
#' @param events A matrix representing the initial events with columns \code{time}, \code{from},
#'   \code{to}, and \code{type} (1 for start, 0 for end, 3 for exogenous changes).
#'   Defaults to an empty 4-column matrix.
#' @param formula_0_1 A one-sided \code{\link[stats]{formula}} specifying the
#'   sufficient statistics for the formation process (\eqn{0 \rightarrow 1}).
#'   The right-hand side must be composed of terms from
#'   \code{\link{redeem_terms}}. For example:
#'   \code{~ inertia() + degree}.
#'   An intercept (\code{~ 1}) is the minimal specification. Defaults to NULL.
#' @param formula_1_0 A one-sided \code{\link[stats]{formula}} specifying the
#'   sufficient statistics for the dissolution process (\eqn{1 \rightarrow 0}).
#'   The right-hand side must be composed of terms from
#'   \code{\link{redeem_terms}}.
#'   An intercept (\code{~ 1}) is the minimal specification. Defaults to NULL.
#' @param coef_0_1 Numeric vector; coefficients for the formation process (\eqn{0 \rightarrow 1}). Defaults to an empty numeric vector.
#' @param coef_1_0 Numeric vector; coefficients for the dissolution process (\eqn{1 \rightarrow 0}). Defaults to an empty numeric vector.
#' @param coef_degree_0_1 Numeric; degree coefficient for the formation process (\eqn{0 \rightarrow 1}). Defaults to 0.
#' @param coef_degree_1_0 Numeric; degree coefficient for the dissolution process (\eqn{1 \rightarrow 0}). Defaults to 0.
#' @param n_events Integer; number of events to simulate. Defaults to 0.
#' @param time Numeric; simulation time limit. Defaults to 0.
#' @param max_events Integer; maximum number of total events. Defaults to 400000.
#' @param n_nodes Integer; the total number of actors in the network.
#' @param verbose Logical; if \code{TRUE}, prints progress information. Defaults to FALSE.
#' @param baseline_0_1 Numeric vector; baseline for the 0 to 1 transition. If the formula for this process
#'   contains an \code{\link{Intercept}} or a \code{\link{degree}} term, then \code{baseline_0_1} should
#'   be a numeric vector with length equal to the number of changepoints, representing
#'   the shifts in the baseline for each interval after the first. If the formula contains
#'   neither, then \code{baseline_0_1} must have length equal to the number of changepoints + 1. Defaults to NULL.
#' @param baseline_1_0 Numeric vector; baseline for the 1 to 0 transition. Similar to \code{baseline_0_1},
#'   its length depends on whether the 1 to 0 formula contains an intercept or degree term. Defaults to NULL.
#' @param simultaneous_interactions Logical; whether to allow simultaneous interactions (i.e. multiple active events for the same actor or dyad at the same time). Defaults to TRUE.
#' @param seed Integer; random seed for simulation. Defaults to 123.
#' @param directed Logical; whether the interaction events are directed. Defaults to FALSE.
#' @note Multi-stream event models are currently not supported in simulation.
#' @return A matrix of simulated events.
#' @examples
#' # Simulate events from a DEM model structure
#' n <- 10
#' f_0_1 <- ~ 1 + inertia(transformation = "identity")
#' f_1_0 <- ~ 1
#'
#' # Simulating events
#' evs <- dem.simulate(
#'   formula_0_1 = f_0_1,
#'   formula_1_0 = f_1_0,
#'   n_nodes = n,
#'   time = 2.0,
#'   coef_0_1 = c(1.0, 0.5),
#'   coef_1_0 = c(-0.5),
#'   seed = 42,
#'   max_events = 100
#' )
#' head(evs)
#' @export
dem.simulate <- function(events = matrix(0, 0, 4),
                         formula_0_1 = NULL,
                         formula_1_0 = NULL,
                         coef_0_1 = numeric(0),
                         coef_1_0 = numeric(0),
                         coef_degree_0_1 = 0,
                         coef_degree_1_0 = 0,
                         n_events = 0,
                         time = 0,
                         max_events = 400000,
                         n_nodes,
                         verbose = FALSE,
                         baseline_0_1 = NULL,
                         baseline_1_0 = NULL,
                         simultaneous_interactions = TRUE,
                         seed = 123,
                         directed = FALSE) {
  set.seed(seed)
  preprocessed <- formula_preprocess(
    formula_1_0 = formula_1_0,
    formula_0_1 = formula_0_1,
    events = events,
    n_nodes = n_nodes,
    directed = directed,
    model_type = "dem",
    simulation = TRUE
  )

  # fixed_effects_0_1 <- if (!is.null(preprocessed$preprocess_0_1)) preprocessed$preprocess_0_1$includes_degrees else FALSE
  # fixed_effects_1_0 <- if (!is.null(preprocessed$preprocess_1_0)) preprocessed$preprocess_1_0$includes_degrees else FALSE
  # fixed_effects <- (fixed_effects_0_1 + fixed_effects_1_0) > 0
  # full_baseline_0_1 <- !fixed_effects_0_1 && length(names(preprocessed$preprocess_0_1$coef_names)) == 0
  # full_baseline_1_0 <- !fixed_effects_1_0 && length(names(preprocessed$preprocess_1_0$coef_names)) == 0

  if (any(!sapply(preprocessed$stream_list, is.null))) {
    stop("Multi-stream simulation is not yet supported.")
  }

  if (!isTRUE(preprocessed$preprocess_0_1$includes_degrees)) {
    if (!is.null(coef_degree_0_1) && (length(coef_degree_0_1) > 1 || any(coef_degree_0_1 != 0))) {
      stop("coef_degree_0_1 provided but degrees are not included in the formula_0_1.")
    }
  }

  if (!isTRUE(preprocessed$preprocess_1_0$includes_degrees)) {
    if (!is.null(coef_degree_1_0) && (length(coef_degree_1_0) > 1 || any(coef_degree_1_0 != 0))) {
      stop("coef_degree_1_0 provided but degrees are not included in the formula_1_0.")
    }
  }

  if (length(baseline_0_1) > 0 && !isTRUE(preprocessed$preprocess_0_1$has_baseline_term)) {
    stop("baseline_0_1 estimates were provided but baseline is not included in the formula_0_1.")
  }

  if (length(baseline_1_0) > 0 && !isTRUE(preprocessed$preprocess_1_0$has_baseline_term)) {
    stop("baseline_1_0 estimates were provided but baseline is not included in the formula_1_0.")
  }
  if ((n_events == 0) && (time == 0)) {
    stop("You need to specify either 'n_events' or 'time'.")
  }
  if ((n_events > 0) && (time > 0)) {
    stop("You need to specify either 'n_events' or 'time', not both.")
  }
  # Match coefficients for both processes
  # If Intercept was removed due to fixed effects, we skip the first unnamed coefficient if provided
  coef_1_0_to_match <- coef_1_0
  if (isTRUE(preprocessed$preprocess_1_0$intercept_removed) && length(coef_1_0) > 0 && (is.null(names(coef_1_0)) || names(coef_1_0)[1] == "")) {
    coef_1_0_to_match <- coef_1_0[-1]
  }
  coef_1_0_augmented <- match_coefficients(
    user_coefs = coef_1_0_to_match,
    internal_names = preprocessed$preprocess_1_0$coef_names,
    internal_keys = names(preprocessed$preprocess_1_0$coef_names)
  )

  coef_0_1_to_match <- coef_0_1
  if (isTRUE(preprocessed$preprocess_0_1$intercept_removed) && length(coef_0_1) > 0 && (is.null(names(coef_0_1)) || names(coef_0_1)[1] == "")) {
    coef_0_1_to_match <- coef_0_1[-1]
  }
  coef_0_1_augmented <- match_coefficients(
    user_coefs = coef_0_1_to_match,
    internal_names = preprocessed$preprocess_0_1$coef_names,
    internal_keys = names(preprocessed$preprocess_0_1$coef_names)
  )


  # Map to full vectors
  coef_0_1_full <- numeric(length(preprocessed$coef_names))
  idx_map_0_1 <- match(preprocessed$preprocess_0_1$coef_names, preprocessed$coef_names)
  coef_0_1_full[idx_map_0_1] <- coef_0_1_augmented

  coef_1_0_full <- numeric(length(preprocessed$coef_names))
  idx_map_1_0 <- match(preprocessed$preprocess_1_0$coef_names, preprocessed$coef_names)
  coef_1_0_full[idx_map_1_0] <- coef_1_0_augmented

  # Detection of changepoints from covariates
  cov_changepoints <- NULL
  if (length(preprocessed$data_list) > 0) {
    for (i in seq_along(preprocessed$data_list)) {
      if (is.list(preprocessed$data_list[[i]])) {
        nms <- names(preprocessed$data_list[[i]])
        if (!is.null(nms)) {
          suppressWarnings(ts <- as.numeric(nms))
          ts <- ts[!is.na(ts) & ts > 1e-10]
          cov_changepoints <- c(cov_changepoints, ts)
        }
      }
    }
  }

  all_changepoints <- sort(unique(c(cov_changepoints, preprocessed$baseline_changepoints_0_1, preprocessed$baseline_changepoints_1_0)))
  if (is.null(all_changepoints)) all_changepoints <- numeric(0)

  if ((n_events == 0) && (time == 0)) {
    if (length(all_changepoints) > 0) {
      time <- max(all_changepoints)
    } else {
      stop("You need to specify either 'n_events' or 'time'.")
    }
  }

  if ((n_events > 0) && (time > 0)) {
    stop("You need to specify either 'n_events' or 'time', not both.")
  }

  window_info <- if (length(preprocessed$window_map) > 0) {
    as.list(stats::setNames(as.numeric(names(preprocessed$window_map)), preprocessed$window_map))
  } else {
    list()
  }

  if (length(all_changepoints) > 0 || !is.null(baseline_0_1) || !is.null(baseline_1_0)) {
    # 1. Align baseline_0_1
    has_intercept_0_1 <- any(tolower(preprocessed$preprocess_0_1$coef_names) == "intercept")
    has_degrees_0_1 <- isTRUE(preprocessed$preprocess_0_1$includes_degrees) || (!is.null(coef_degree_0_1) && any(coef_degree_0_1 != 0))
    full_baseline_0_1 <- !(has_intercept_0_1 || has_degrees_0_1)

    ref_tc_0_1 <- if (!is.null(preprocessed$baseline_changepoints_0_1)) {
      preprocessed$baseline_changepoints_0_1
    } else {
      numeric(0)
    }

    if (!is.null(baseline_0_1)) {
      expected_len_0_1 <- if (full_baseline_0_1) length(ref_tc_0_1) + 1 else length(ref_tc_0_1)
      if (full_baseline_0_1 && length(baseline_0_1) == length(ref_tc_0_1)) {
        baseline_0_1 <- c(0, baseline_0_1)
        expected_len_0_1 <- length(ref_tc_0_1) + 1
      }
      if (length(baseline_0_1) != expected_len_0_1) stop("Incorrect length for baseline_0_1")
      baseline_0_1_full <- if (full_baseline_0_1) baseline_0_1 else c(0, baseline_0_1)
      full_tc <- c(-1e-10, ref_tc_0_1, Inf)
      slice_indices <- findInterval(c(0, all_changepoints), full_tc)
      baseline_0_1 <- baseline_0_1_full[slice_indices]
    } else {
      baseline_0_1 <- rep(0, length(all_changepoints) + 1)
    }

    # 2. Align baseline_1_0
    has_intercept_1_0 <- any(tolower(preprocessed$preprocess_1_0$coef_names) == "intercept")
    has_degrees_1_0 <- isTRUE(preprocessed$preprocess_1_0$includes_degrees) || (!is.null(coef_degree_1_0) && any(coef_degree_1_0 != 0))
    full_baseline_1_0 <- !(has_intercept_1_0 || has_degrees_1_0)

    ref_tc_1_0 <- if (!is.null(preprocessed$baseline_changepoints_1_0)) {
      preprocessed$baseline_changepoints_1_0
    } else {
      numeric(0)
    }

    if (!is.null(baseline_1_0)) {
      expected_len_1_0 <- if (full_baseline_1_0) length(ref_tc_1_0) + 1 else length(ref_tc_1_0)
      if (full_baseline_1_0 && length(baseline_1_0) == length(ref_tc_1_0)) {
        baseline_1_0 <- c(0, baseline_1_0)
        expected_len_1_0 <- length(ref_tc_1_0) + 1
      }
      if (length(baseline_1_0) != expected_len_1_0) stop("Incorrect length for baseline_1_0")
      baseline_1_0_full <- if (full_baseline_1_0) baseline_1_0 else c(0, baseline_1_0)
      full_tc <- c(-1e-10, ref_tc_1_0, Inf)
      slice_indices <- findInterval(c(0, all_changepoints), full_tc)
      baseline_1_0 <- baseline_1_0_full[slice_indices]
    } else {
      baseline_1_0 <- rep(0, length(all_changepoints) + 1)
    }

    if (verbose) {
      print_simulation_info(
        preprocessed = preprocessed,
        coef_0_1_full = coef_0_1_full,
        coef_1_0_full = coef_1_0_full,
        baseline_0_1 = baseline_0_1,
        baseline_1_0 = baseline_1_0,
        all_changepoints = all_changepoints,
        time = time
      )
    }

    res <- simulate_from_empty_timevarying(
      terms = preprocessed$term_names,
      data_list = preprocessed$data_list,
      transformations = preprocessed$transformation_list,
      n_nodes = preprocessed$n_nodes,
      verbose = verbose,
      directed = directed,
      coef_0_1 = coef_0_1_full,
      coef_1_0 = coef_1_0_full,
      degree_coef_0_1 = coef_degree_0_1,
      degree_coef_1_0 = coef_degree_1_0,
      time_changepoints = all_changepoints,
      baseline_0_1 = baseline_0_1,
      baseline_1_0 = baseline_1_0,
      simultaneous_interactions = simultaneous_interactions,
      n_events = n_events,
      max_events = max_events,
      time = time,
      seed = seed,
      window_info = window_info
    )
    data_tmp <- res[length(res)]
    res <- res[-length(res)]

    res <- do.call(rbind, res)
    return(res)
  } else {
    if (verbose) {
      print_simulation_info(
        preprocessed = preprocessed,
        coef_0_1_full = coef_0_1_full,
        coef_1_0_full = coef_1_0_full,
        baseline_0_1 = 0,
        baseline_1_0 = 0,
        all_changepoints = numeric(0)
      )
    }
    return(
      simulate_from_empty(
        terms = preprocessed$term_names,
        data_list = preprocessed$data_list,
        transformations = preprocessed$transformation_list,
        n_nodes = preprocessed$n_nodes,
        verbose = verbose,
        directed = directed,
        coef_0_1 = coef_0_1_full,
        coef_1_0 = coef_1_0_full,
        degree_coef_0_1 = coef_degree_0_1,
        degree_coef_1_0 = coef_degree_1_0,
        n_events = n_events,
        time = time,
        max_events = max_events,
        simultaneous_interactions = simultaneous_interactions,
        seed = seed,
        window_info = window_info
      )
    )
  }
}
