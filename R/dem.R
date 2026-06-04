#' Durational Event Model (DEM) Estimation
#'
#' Estimates a Durational Event Model (DEM) for relational event sequences
#' where interactions have a duration.
#'
#' @details
#' The Durational Event Model (DEM) is a general framework for analyzing durational
#' events, extending standard Relational Event Models (REM) by decoupling the
#' modeling of event incidence from event duration. It characterizes the dynamics
#' via two separate continuous-time counting processes:
#'
#' \describe{
#'   \item{Formation Process (\eqn{0 \rightarrow 1})}{Counts the number of times that
#'   actor pair \eqn{(i,j)} starts an interaction up to time \eqn{t}. The incidence intensity
#'   is denoted by \eqn{\lambda_{i,j}^{0\rightarrow 1}(t | \mathscr{H}_t)}.}
#'   \item{Dissolution Process (\eqn{1 \rightarrow 0})}{Counts the number of times that
#'   actor pair \eqn{(i,j)} stops interacting up to time \eqn{t}. The dissolution intensity
#'   is denoted by \eqn{\lambda_{i,j}^{1\rightarrow 0}(t | \mathscr{H}_t)}.}
#' }
#'
#' Under the assumption that the processes are non-homogeneous Poisson processes, the intensities
#' are modeled as:
#' \deqn{\lambda_{i,j}^{0\rightarrow 1}(t | \mathscr{H}_t, \theta^{0\rightarrow 1}) = \exp(s_{i,j}^{0\rightarrow 1}(\mathscr{H}_t)^\top \alpha^{0\rightarrow 1} + \beta_i^{0\rightarrow 1} + \beta_j^{0\rightarrow 1} + f(t, \gamma^{0\rightarrow 1}))}
#' \deqn{\lambda_{i,j}^{1\rightarrow 0}(t | \mathscr{H}_t, \theta^{1\rightarrow 0}) = \exp(s_{i,j}^{1\rightarrow 0}(\mathscr{H}_t)^\top \alpha^{1\rightarrow 0} + \beta_i^{1\rightarrow 0} + \beta_j^{1\rightarrow 0} + f(t, \gamma^{1\rightarrow 0}))}
#' where:
#' \itemize{
#'   \item \eqn{s_{i,j}(\mathscr{H}_t)} is a vector of dynamic network statistics
#'     capturing the history of past interactions \eqn{\mathscr{H}_t}.
#'   \item \eqn{\alpha} is a parameter vector determining the covariate effects.
#'   \item \eqn{\beta_i} and \eqn{\beta_j} are actor-specific sociality/popularity
#'     parameters (degree correction) capturing actor heterogeneity.
#'   \item \eqn{f(t, \gamma)} is a piecewise-constant step function modeling
#'     temporal baseline fluctuations across a set of changepoints.
#' }
#'
#' To satisfy the Feller criterion and ensure that the continuous-time
#' counting process remains non-explosive, count-based network statistics
#' (such as inertia or common partners) are typically log-transformed on
#' the \eqn{\log(x + 1)} scale.
#'
#' @section Scalable Estimation Algorithm:
#' The likelihood of the model is separable with respect to
#' \eqn{\theta^{0\rightarrow 1}} and \eqn{\theta^{1\rightarrow 0}}, allowing
#' independent estimation of the incidence and duration components.
#' Traditional maximum likelihood estimation via standard Newton-Raphson
#' requires computing and inverting an \eqn{O(N^2)} Hessian matrix, which
#' is computationally prohibitive for larger networks. To bypass this, the
#' \code{redeem} package implements a highly scalable block-coordinate
#' ascent algorithm that separates parameter updates:
#' \enumerate{
#'   \item \strong{Step 1}: Update covariate parameters \eqn{\alpha} using
#'     a standard Newton-Raphson update.
#'   \item \strong{Step 2}: Update high-dimensional actor popularity
#'     baselines \eqn{\beta} using Minorization-Maximization
#'     (MM) steps, avoiding explicit matrix inversion.
#'   \item \strong{Step 3}: Update baseline step function parameters
#'     \eqn{\gamma} via a closed-form step.
#' }
#' More information is provided in Fritz et at. (2026).
#'
#' @section Semiparametric Baseline:
#' When \code{semiparametric = TRUE}, the baseline rates for both the formation
#' (\eqn{0 \rightarrow 1}) and dissolution (\eqn{1 \rightarrow 0}) processes are
#' left completely unspecified. Both processes are estimated as separate Cox proportional
#' hazards models using the \code{survival} package. In this path:
#' \itemize{
#'   \item For the formation process, each start event (1) is treated as a failure, and
#'     all inactive dyads at that time constitute the risk set.
#'   \item For the dissolution process, each end event (0) is treated as a failure, and
#'     all currently active interactions constitute the risk set.
#'   \item The exact waiting times (durations of the active and inactive states) are
#'     conditioned out, and estimation is based solely on the ordering of events and
#'     the relative dyadic intensities at each transition time.
#'   \item This approach is highly robust to arbitrary temporal fluctuations in baseline rates
#'     since no piecewise-constant temporal baselines or changepoints need to be specified.
#'   \item \strong{Limitations}: This path does \strong{not} support the specialized
#'     scalable estimation of sender/receiver popularity effects (\code{degree}) or
#'     piecewise-constant temporal baselines.
#' }
#'
#' @references
#' Fritz, C., Rastelli, R., Fop, M., & Caimo, A. (2026). Scalable Durational Event Models:
#' Application to Physical and Digital Interactions. arXiv:2504.00049.
#'
#' @param events A matrix of events with columns \code{time}, \code{from}, \code{to},
#'   and \code{type} (1 for start, 0 for end, 3 for exogenous changes).
#' @param training_start Numeric; the time point at which to start the estimation.
#'   Defaults to 0.
#' @param exogenous_end Numeric; the exogenous time point at which the observational period ends. 
#'   Defaults to NULL, which implies that time when the final event was observed is taken as the end of the observational period.
#'
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
#' @param n_nodes Integer; the total number of actors in the network. If \code{NULL} (default), it is automatically identified based on the actors in the \code{events} set.
#' @param directed Logical; whether the interaction events are directed.
#'   Defaults to FALSE.
#' @param estimate_0_1 Logical; whether to estimate the formation process.
#'   Defaults to NULL, in which case it is estimated if \code{formula_0_1} is provided.
#' @param estimate_1_0 Logical; whether to estimate the dissolution process.
#'   Defaults to NULL, in which case it is estimated if \code{formula_1_0} is provided.
#' @param coef_0_1 Numeric vector; initial coefficients for the formation model.
#'   If provided, this must be a concatenated vector of:
#'   \enumerate{
#'     \item Core coefficients: values for sufficient statistics in the formula.
#'     \item Degree coefficients (if \code{degree} is in the formula): a vector of length \code{n_nodes} (undirected) or \code{2 * n_nodes} (directed, sender effects first then receiver effects).
#'     \item Baseline coefficients (if temporal changepoints are present): a vector of length equal to the number of baseline intervals (equal to number of changepoints if an intercept/degree is present, or changepoints + 1 if neither is present).
#'   }
#'   Defaults to NULL, in which case default starting values are automatically computed.
#' @param coef_1_0 Numeric vector; initial coefficients for the dissolution model.
#'   If provided, this must be a concatenated vector of:
#'   \enumerate{
#'     \item Core coefficients: values for sufficient statistics in the formula.
#'     \item Degree coefficients (if \code{degree} is in the formula): a vector of length \code{n_nodes} (undirected) or \code{2 * n_nodes} (directed, sender effects first then receiver effects).
#'     \item Baseline coefficients (if temporal changepoints are present): a vector of length equal to the number of baseline intervals (equal to number of changepoints if an intercept/degree is present, or changepoints + 1 if neither is present).
#'   }
#'   Defaults to NULL, in which case default starting values are automatically computed.
#' @param semiparametric Logical; whether to use a semiparametric baseline.
#'   Defaults to FALSE. See the 'Semiparametric Baseline' section for details.
#' @param simultaneous_interactions Logical; whether to allow simultaneous interactions
#'   (i.e. multiple active events for the same actor or dyad at the same time).
#'   Defaults to TRUE.
#' @param control A list of control parameters from \code{\link{control.redeem}}.
#'   Defaults to \code{control.redeem()}.
#'
#' @return An object of class \code{\link{dem_object}} containing model estimates,
#'   log-likelihoods, and preprocessed data. See \code{\link{dem_object}}
#'   for details on the components of the returned object and S3 methods.
#'
#' @examples
#' # Simulate some durational data
#' n <- 20
#' events <- matrix(c(
#'   1.2, 1, 5, 1,
#'   2.5, 1, 5, 0,
#'   3.1, 2, 8, 1,
#'   4.4, 2, 8, 0
#' ), ncol = 4, byrow = TRUE)
#' colnames(events) <- c("time", "from", "to", "type")
#'
#' # Estimate a simple DEM
#' fit <- dem(
#'   events = events,
#'   n_nodes = n,
#'   formula_0_1 = ~1,
#'   formula_1_0 = ~1,
#'   control = control.redeem(estimate = "Blockwise")
#' )
#' summary(fit)
#' @export
dem <- function(events,
                training_start = 0,
                exogenous_end = NULL,
                formula_0_1 = NULL,
                formula_1_0 = NULL,
                n_nodes = NULL,
                directed = FALSE,
                estimate_0_1 = NULL,
                estimate_1_0 = NULL,
                coef_0_1 = NULL,
                coef_1_0 = NULL,
                semiparametric = FALSE,
                simultaneous_interactions = TRUE,
                control = control.redeem()) {
  call <- match.call()

  processed_actors <- process_event_actors(events, n_nodes, directed = directed)
  events <- processed_actors$events
  n_nodes <- processed_actors$n_nodes

  # Convert to matrix if it's a data.frame or data.table to ensure consistent indexing
  if (!is.matrix(events)) events <- as.matrix(events)

  if (ncol(events) < 4) {
    stop(sprintf("Event matrix must have at least 4 columns (time, from, to, type). Found only %d columns.", ncol(events)))
  }

  if (!is.logical(simultaneous_interactions)) {
    stop("simultaneous_interactions must be logical.")
  }

  if (!is.null(exogenous_end)) {
    if (nrow(events) > 0 && max(events[, 1]) > exogenous_end) {
      warning("exogenous_end is before the last event. Truncating events.")
      events <- events[events[, 1] <= exogenous_end, , drop = FALSE]
    }
  }

  # Detect columns by name if possible, otherwise assume standard order (1=time, 2=from, 3=to, 4=type)
  b_time_val <- if (!is.null(control$build_time)) control$build_time else 0
  effective_start <- min(training_start, b_time_val, events[, 1], na.rm = TRUE)

  if (control$check_matrix) {
    res_check <- check_matrix(events, return_matrix = TRUE, start_time = effective_start)
    if (!is.null(attr(res_check, "has_errors")) && attr(res_check, "has_errors")) {
      stop("The provided event data contains fatal errors (e.g. overlaps). Estimation aborted. Please correct the issues listed in the warnings above.")
    }
    if (is.matrix(res_check)) {
      events <- res_check
    }
  }

  # Extract control parameters
  it_max <- rep(control$it_max, length.out = 2)
  tol <- rep(control$tol, length.out = 2)
  accelerated <- rep(control$accelerated, length.out = 2)
  verbose <- control$verbose
  weighting <- control$weighting
  subsample <- control$subsample
  return_data <- control$return_data
  save_hist <- control$save_hist
  estimate <- control$estimate
  build_time <- control$build_time
  use_glm <- control$use_glm
  legacy <- control$legacy
  inf_unidentifiable <- control$inf_unidentifiable

  if (ncol(events) != 4) {
    stop("Events matrix for dem() must have 4 columns: time, from, to, and type (1 for start, 0 for end).")
  }
  colnames(events) <- c("time", "from", "to", "type")

  # Handle exogenous observation end early to ensure it's captured in event history and formula preprocessing
  if (!is.null(exogenous_end)) {
    events <- rbind(events, matrix(c(exogenous_end, 1, 2, 3), nrow = 1))
    events <- events[order(events[, 1]), ]
  }

  number_0_1 <- sum(events[, 4] == 1)
  number_1_0 <- sum(events[, 4] == 0)
  event_numbers <- c(number_0_1, number_1_0)

  # 1. Preprocess formulas
  preprocessed <- formula_preprocess(
    formula_0_1 = formula_0_1,
    formula_1_0 = formula_1_0,
    model_type = "dem",
    events = events,
    n_nodes = n_nodes,
    directed = directed
  )

  if (is.null(estimate_0_1)) {
    estimate_0_1 <- !is.null(formula_0_1)
  }
  if (is.null(estimate_1_0)) {
    estimate_1_0 <- !is.null(formula_1_0)
  }

  if (is.null(formula_0_1)) {
    estimate_0_1 <- FALSE
  }
  if (is.null(formula_1_0)) {
    estimate_1_0 <- FALSE
  }


  # 2. Extract shared information
  time_changepoints_0_1 <- preprocessed$baseline_changepoints_0_1
  time_changepoints_1_0 <- preprocessed$baseline_changepoints_1_0
  all_changepoints <- sort(unique(c(time_changepoints_0_1, time_changepoints_1_0)))

  # Ensure study period covers all changepoints
  if (is.null(exogenous_end)) {
    last_event_time <- if (nrow(events) > 0) max(events[, 1], na.rm = TRUE) else 0
    if (length(all_changepoints) > 0) {
      max_cp <- max(all_changepoints)
      if (max_cp > last_event_time) {
        warning(paste("Study end (max event time =", last_event_time, ") is before the last baseline changepoint (", max_cp, "). Extending study period to", max_cp, ". For durational models, consider providing 'exogenous_end' explicitly."))
        exogenous_end <- max_cp
      } else {
        exogenous_end <- last_event_time
      }
    } else {
      exogenous_end <- last_event_time
    }
  }

  # If exogenous changes to the baseline intensity are provided, we add them as artificial events with type 3
  if (length(all_changepoints) > 0) {
    if (estimate == "NR") {
      events <- rbind(events, cbind(all_changepoints, 1, 2, 3))
      events <- events[order(events[, 1]), ]
      preprocessed$events <- events
    }
  }

  time_labels_0_1 <- preprocessed$baseline_labels_0_1
  time_labels_1_0 <- preprocessed$baseline_labels_1_0

  if (training_start != 0) {
    if (nrow(events) > 0 && training_start >= max(events[, 1])) {
      stop("The proposed starting time is after the last event.")
    }
    # If we start not from the beginning all change points that are before
    # the training will be deleted
    include_0_1 <- time_changepoints_0_1 > training_start
    time_changepoints_0_1 <- time_changepoints_0_1[include_0_1]
    if (!is.null(time_labels_0_1)) time_labels_0_1 <- time_labels_0_1[include_0_1]

    include_1_0 <- time_changepoints_1_0 > training_start
    time_changepoints_1_0 <- time_changepoints_1_0[include_1_0]
    if (!is.null(time_labels_1_0)) time_labels_1_0 <- time_labels_1_0[include_1_0]
  }

  labels_changepoints_0_1 <- if (!is.null(time_labels_0_1)) time_labels_0_1 else as.character(time_changepoints_0_1)
  labels_changepoints_1_0 <- if (!is.null(time_labels_1_0)) time_labels_1_0 else as.character(time_changepoints_1_0)
  if (!(estimate != "Blockwise") && weighting) {
    weighting <- FALSE
    if (verbose) {
      message(
        paste(
          "Weighting is currently only supported for the current model when it is estimated with Newton Raphon Methods.\n",
          "Therefore, we have set the weighting parameter to FALSE."
        )
      )
    }
  }

  fixed_effects_0_1 <- if (!is.null(preprocessed$preprocess_0_1)) preprocessed$preprocess_0_1$includes_degrees else FALSE
  fixed_effects_1_0 <- if (!is.null(preprocessed$preprocess_1_0)) preprocessed$preprocess_1_0$includes_degrees else FALSE
  fixed_effects <- (fixed_effects_0_1 + fixed_effects_1_0) > 0

  if (!fixed_effects && is.null(time_changepoints_0_1) && is.null(time_changepoints_1_0) && estimate == "Blockwise") {
    estimate <- "NR"
  }
  # If the baseline intensity is a constant and no fixed effects are needed,
  # the estimation can be done with glm
  if (verbose) {
    cat("Calling preprocess with:", "\n")
    cat("edgelist dim: ", paste(dim(as.matrix(preprocessed$events)), collapse = "x"), "\n")
    cat("terms: ", paste(preprocessed$term_names, collapse = ", "), "\n")
    cat("n_nodes: ", n_nodes, "\n")
    cat("simultaneous_interactions: ", simultaneous_interactions, "\n")
  }

  if (is.null(build_time)) {
    build_time <- 0
  }
  build_time <- max(training_start, build_time)

  max_time <- if (nrow(events) > 0) max(events[, 1]) else 0

  preprocessed_tmp <- preprocess_multi_stream(
    preprocessed = preprocessed,
    n_nodes = n_nodes,
    verbose = verbose,
    directed = directed,
    simultaneous_interactions = simultaneous_interactions,
    build_time = build_time,
    max_time = max_time,
    model_type = "dem"
  )

  if (!is.null(build_time) && build_time > 0) {
    preprocessed_tmp <- preprocessed_tmp[preprocessed_tmp$time >= build_time, ]
  }

  if (is.null(preprocessed_tmp) || ncol(preprocessed_tmp) == 0 || nrow(preprocessed_tmp) == 0) {
    stop("Preprocessing failed to produce valid data intervals. Check if your event sequence is valid (e.g., each start has an end).")
  }
  # Ensure we use the actual coefficient names, not the names attribute which might be NULL
  coef_names_vec <- preprocessed$coef_names

  expected_names <- c(
    "time_new",
    "time",
    "pair_id",
    "status",
    "event",
    "from",
    "to",
    "from_avail",
    "to_avail",
    coef_names_vec
  )

  if (ncol(preprocessed_tmp) != length(expected_names)) {
    warning(paste("Dimension mismatch between preprocessed data (", ncol(preprocessed_tmp), ") and coefficient names (", length(expected_names), "). Adjusting names."))

    if (ncol(preprocessed_tmp) <= length(expected_names)) {
      colnames(preprocessed_tmp) <- expected_names[seq_len(ncol(preprocessed_tmp))]
    } else {
      colnames(preprocessed_tmp) <- c(expected_names, paste0("V", (length(expected_names) + 1):ncol(preprocessed_tmp)))
    }
  } else {
    colnames(preprocessed_tmp) <- expected_names
  }

  preprocessed_tmp <- data.table::as.data.table(preprocessed_tmp)
  preprocessed_tmp[, `:=`(
    diff = time_new - time
  )]
  preprocessed_tmp[diff <= 0, diff := NA]
  preprocessed_tmp[, `:=`(
    offset = log(diff),
    avail = (from_avail + to_avail) == 2
  )]
  preprocessed_tmp[is.infinite(offset), event := NA]

  matching <- c("event", "offset", "status", names(preprocessed$coef_names))
  if (fixed_effects) {
    matching <- c(matching, "from", "to")
  }
  preprocessed_tmp <- preprocessed_tmp[!is.na(event) & !is.na(from) & !is.na(to), ]

  formula_1_0_names <- names(preprocessed$preprocess_1_0$coef_names)
  # If Intercept was removed due to degrees and no other core terms exist, use ~1 but it will be ignored if missing from data
  formula_1_0_new <- if (length(formula_1_0_names) > 0) {
    formula(paste("~", paste(collapse = "+ ", formula_1_0_names)))
  } else {
    ~1
  }

  formula_0_1_names <- names(preprocessed$preprocess_0_1$coef_names)
  formula_0_1_new <- if (length(formula_0_1_names) > 0) {
    formula(paste("~", paste(collapse = "+ ", formula_0_1_names)))
  } else {
    ~1
  }


  if (!simultaneous_interactions) {
    matching <- c(matching, "avail")
  }

  if (estimate != "Blockwise") {
    if (!is.null(time_changepoints_0_1) && length(time_changepoints_0_1) > 0) {
      preprocessed_tmp$time_cat_0_1 <- cut(preprocessed_tmp$time_new, breaks = c(-1, time_changepoints_0_1, Inf), labels = c("Beg", labels_changepoints_0_1))
      matching <- c(matching, "time_cat_0_1")
    }
    if (!is.null(time_changepoints_1_0) && length(time_changepoints_1_0) > 0) {
      preprocessed_tmp$time_cat_1_0 <- cut(preprocessed_tmp$time_new, breaks = c(-1, time_changepoints_1_0, Inf), labels = c("Beg", labels_changepoints_1_0))
      matching <- c(matching, "time_cat_1_0")
    }
  }

  matching <- c(matching, "pair_id")
  preprocessed_tmp$weight <- 1
  data.table::setDT(preprocessed_tmp)
  if (weighting) {
    preprocessed_tmp <- preprocessed_tmp[, .(weight = sum(weight)), by = c(matching)]
  }
  # Prepare transition-specific data
  data_1_0 <- preprocessed_tmp[status == 1]
  data_0_1 <- preprocessed_tmp[status == 0]
  if (!simultaneous_interactions) {
    data_0_1 <- data_0_1[avail == TRUE]
  }
  rm(preprocessed_tmp)

  if (estimate != "Blockwise") {
    if (fixed_effects) {
      if (directed) {
        sender_names <- paste0("sender_", seq_len(n_nodes))
        receiver_names <- paste0("receiver_", seq_len(n_nodes))
        fixed_effect_names <- c(sender_names, receiver_names)

        # Transition 0-1
        fixed_effect_0_1_s <- matrix(0, nrow = nrow(data_0_1), ncol = n_nodes)
        fixed_effect_0_1_s[cbind(seq_len(nrow(data_0_1)), data_0_1$from)] <- 1
        colnames(fixed_effect_0_1_s) <- sender_names

        fixed_effect_0_1_r <- matrix(0, nrow = nrow(data_0_1), ncol = n_nodes)
        fixed_effect_0_1_r[cbind(seq_len(nrow(data_0_1)), data_0_1$to)] <- 1
        colnames(fixed_effect_0_1_r) <- receiver_names

        data_0_1 <- cbind(data_0_1, fixed_effect_0_1_s, fixed_effect_0_1_r)
        formula_0_1_new <- update(formula_0_1_new, new = paste("~ . -Intercept +", paste(fixed_effect_names, collapse = "+")))

        # Transition 1-0
        fixed_effect_1_0_s <- matrix(0, nrow = nrow(data_1_0), ncol = n_nodes)
        fixed_effect_1_0_s[cbind(seq_len(nrow(data_1_0)), data_1_0$from)] <- 1
        colnames(fixed_effect_1_0_s) <- sender_names

        fixed_effect_1_0_r <- matrix(0, nrow = nrow(data_1_0), ncol = n_nodes)
        fixed_effect_1_0_r[cbind(seq_len(nrow(data_1_0)), data_1_0$to)] <- 1
        colnames(fixed_effect_1_0_r) <- receiver_names

        data_1_0 <- cbind(data_1_0, fixed_effect_1_0_s, fixed_effect_1_0_r)
        formula_1_0_new <- update(formula_1_0_new, new = paste("~ . -Intercept +", paste(fixed_effect_names, collapse = "+")))
      } else {
        fixed_effect_names <- paste0("effect_", seq_len(n_nodes))

        fixed_effect_0_1 <- matrix(0, nrow = nrow(data_0_1), ncol = n_nodes)
        fixed_effect_0_1[cbind(rep(seq_len(nrow(fixed_effect_0_1)), times = 2), c(data_0_1$from, data_0_1$to))] <- 1
        colnames(fixed_effect_0_1) <- fixed_effect_names
        data_0_1 <- cbind(data_0_1, fixed_effect_0_1)
        formula_0_1_new <- update(formula_0_1_new, new = paste("~ . -Intercept +", paste(fixed_effect_names, collapse = "+")))

        fixed_effect_1_0 <- matrix(0, nrow = nrow(data_1_0), ncol = n_nodes)
        fixed_effect_1_0[cbind(rep(seq_len(nrow(fixed_effect_1_0)), times = 2), c(data_1_0$from, data_1_0$to))] <- 1
        colnames(fixed_effect_1_0) <- fixed_effect_names
        data_1_0 <- cbind(data_1_0, fixed_effect_1_0)
        formula_1_0_new <- update(formula_1_0_new, new = paste("~ . -Intercept +", paste(fixed_effect_names, collapse = "+")))
      }
    }
    if (!is.null(time_changepoints_0_1) && length(time_changepoints_0_1) > 0) {
      if (!"time_cat" %in% names(data_0_1)) {
        data_0_1$time_cat <- if ("time_cat_0_1" %in% names(data_0_1)) data_0_1$time_cat_0_1 else cut(data_0_1$time_new, breaks = c(-1, time_changepoints_0_1, Inf), labels = c("Beg", labels_changepoints_0_1))
      }
      data_0_1 <- cbind(data_0_1, model.matrix(~ -1 + time_cat, data_0_1))

      full_baseline_0_1 <- !fixed_effects && length(names(preprocessed$preprocess_0_1$coef_names)) == 0
      if (full_baseline_0_1) {
        time_cat_names <- paste0("time_cat", levels(data_0_1$time_cat))
        formula_0_1_new <- formula(paste("~ -1 +", paste(time_cat_names, collapse = " + ")))
      } else {
        time_cat_names <- paste0("time_cat", levels(data_0_1$time_cat)[-1])
        if (length(time_cat_names) > 0) {
          formula_0_1_new <- update(formula_0_1_new, new = paste("~ . +", paste(time_cat_names, collapse = "+")))
        }
      }
    }
    if (!is.null(time_changepoints_1_0) && length(time_changepoints_1_0) > 0) {
      if (!"time_cat" %in% names(data_1_0)) {
        data_1_0$time_cat <- if ("time_cat_1_0" %in% names(data_1_0)) data_1_0$time_cat_1_0 else cut(data_1_0$time_new, breaks = c(-1, time_changepoints_1_0, Inf), labels = c("Beg", labels_changepoints_1_0))
      }
      data_1_0 <- cbind(data_1_0, model.matrix(~ -1 + time_cat, data_1_0))

      full_baseline_1_0 <- !fixed_effects && length(names(preprocessed$preprocess_1_0$coef_names)) == 0
      if (full_baseline_1_0) {
        time_cat_names <- paste0("time_cat", levels(data_1_0$time_cat))
        formula_1_0_new <- formula(paste("~ -1 +", paste(time_cat_names, collapse = " + ")))
      } else {
        time_cat_names <- paste0("time_cat", levels(data_1_0$time_cat)[-1])
        if (length(time_cat_names) > 0) {
          formula_1_0_new <- update(formula_1_0_new, new = paste("~ . +", paste(time_cat_names, collapse = "+")))
        }
      }
    }
    indicators_0_1 <- indicators_1_0 <- NULL
  } else {
    indicators_0_1 <- match(attr(terms(formula_0_1_new), "term.labels"), names(data_0_1))
    indicators_1_0 <- match(attr(terms(formula_1_0_new), "term.labels"), names(data_1_0))
  }

  # Transition 0-1 (Start)
  start_time <- Sys.time()
  if (estimate_0_1) {
    model_0_1 <- estimate_transition(
      data = data_0_1,
      formula_original = formula_0_1,
      formula_new = formula_0_1_new,
      indicators = indicators_0_1,
      n_nodes = n_nodes,
      estimate_method = estimate,
      it_max = it_max[1],
      tol = tol[1],
      accelerated = accelerated[1],
      subsample = subsample,
      verbose = verbose,
      estimate_degree = fixed_effects,
      directed = directed,
      semiparametric = semiparametric,
      labels_changepoints = labels_changepoints_0_1,
      time_changepoints = time_changepoints_0_1,
      coef_init = coef_0_1,
      model_type = "dem",
      process = "0-1",
      return_data = return_data,
      save_hist = save_hist,
      use_glm = use_glm,
      legacy = legacy,
      inf_unidentifiable = inf_unidentifiable,
      events = events
    )
  } else {
    model_0_1 <- NULL
  }

  # Transition 1-0 (End)
  if (estimate_1_0) {
    model_1_0 <- estimate_transition(
      data = data_1_0,
      formula_original = formula_1_0,
      formula_new = formula_1_0_new,
      indicators = indicators_1_0,
      n_nodes = n_nodes,
      estimate_method = estimate,
      it_max = it_max[2],
      tol = tol[2],
      accelerated = accelerated[2],
      subsample = subsample,
      verbose = verbose,
      estimate_degree = fixed_effects,
      directed = directed,
      semiparametric = semiparametric,
      labels_changepoints = labels_changepoints_1_0,
      time_changepoints = time_changepoints_1_0,
      coef_init = coef_1_0,
      model_type = "dem",
      process = "1-0",
      return_data = return_data,
      save_hist = save_hist,
      use_glm = use_glm,
      legacy = legacy,
      inf_unidentifiable = inf_unidentifiable,
      events = events
    )
  } else {
    model_1_0 <- NULL
  }
  end_time <- Sys.time()

  res <- list(
    call = call,
    event_numbers = event_numbers,
    model_1_0 = model_1_0,
    model_0_1 = model_0_1,
    events = events,
    formula_0_1 = formula_0_1,
    formula_1_0 = formula_1_0,
    n_nodes = n_nodes,
    simultaneous_interactions = simultaneous_interactions,
    directed = directed,
    training_start = training_start,
    build_time = build_time,
    max_time = max_time,
    exogenous_end = exogenous_end,
    time_changepoints = preprocessed$baseline_changepoints_0_1,
    labels_changepoints = preprocessed$baseline_labels_0_1,
    subsample = subsample,
    return_data = return_data,
    runtime = end_time - start_time,
    window_map = preprocessed$window_map,
    preprocessed = preprocessed
  )
  # Apply descriptive names to results
  if (!is.null(res$model_0_1$coefficients)) {
    mapping_0_1 <- preprocessed$preprocess_0_1$coef_names
    curr_names <- names(res$model_0_1$coefficients)
    if (!is.null(curr_names)) {
      new_names <- curr_names
      for (i in seq_along(mapping_0_1)) {
        new_names[new_names == names(mapping_0_1)[i]] <- mapping_0_1[i]
      }
      names(res$model_0_1$coefficients) <- new_names

      # Cascade renaming to other components
      n_core <- length(res$model_0_1$est_core)
      if (n_core > 0) names(res$model_0_1$est_core) <- new_names[seq_len(n_core)]

      if (!is.null(res$model_0_1$covariance)) {
        # Use names that match the actual dimensions of the covariance matrix
        # This handles cases where terms were dropped due to singularity
        actual_names <- colnames(res$model_0_1$covariance)
        if (!is.null(actual_names)) {
          descriptive_cov_names <- actual_names
          for (i in seq_along(mapping_0_1)) {
            descriptive_cov_names[descriptive_cov_names == names(mapping_0_1)[i]] <- mapping_0_1[i]
          }
          colnames(res$model_0_1$covariance) <- rownames(res$model_0_1$covariance) <- descriptive_cov_names
        }
      }
    }
  }
  if (!is.null(res$model_1_0$coefficients)) {
    mapping_1_0 <- preprocessed$preprocess_1_0$coef_names
    curr_names <- names(res$model_1_0$coefficients)
    if (!is.null(curr_names)) {
      new_names <- curr_names
      for (i in seq_along(mapping_1_0)) {
        new_names[new_names == names(mapping_1_0)[i]] <- mapping_1_0[i]
      }
      names(res$model_1_0$coefficients) <- new_names

      # Cascade renaming to other components
      n_core <- length(res$model_1_0$est_core)
      if (n_core > 0) names(res$model_1_0$est_core) <- new_names[seq_len(n_core)]

      if (!is.null(res$model_1_0$covariance)) {
        actual_names <- colnames(res$model_1_0$covariance)
        if (!is.null(actual_names)) {
          descriptive_cov_names <- actual_names
          for (i in seq_along(mapping_1_0)) {
            descriptive_cov_names[descriptive_cov_names == names(mapping_1_0)[i]] <- mapping_1_0[i]
          }
          colnames(res$model_1_0$covariance) <- rownames(res$model_1_0$covariance) <- descriptive_cov_names
        }
      }
    }
  }

  class(res) <- "dem"
  return(res)
}
