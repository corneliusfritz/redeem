#' Relational Event Model (REM) Estimation
#'
#' Estimates a Relational Event Model (REM) for network data, focusing on the
#' incidence of discrete events between pairs of actors.
#' See \code{\link{dem}} for the full Durational Event Model, which extends
#' the REM to handle interactions with non-negligible duration.
#'
#' @section Model Formulation:
#' The Relational Event Model characterizes the instantaneous rate at which
#' actor pair \eqn{(i,j)} initiates an event. Under the log-linear
#' specification, the event intensity at time \eqn{t} is:
#' \deqn{\lambda_{i,j}(t \mid \mathscr{H}_t,\, \theta) =
#'   \exp\!\bigl(s_{i,j}(\mathscr{H}_t)^\top \alpha +
#'   \beta_i + \beta_j + f(t, \gamma)\bigr)}
#' where:
#' \itemize{
#'   \item \eqn{s_{i,j}(\mathscr{H}_t)} is a vector of sufficient statistics
#'     computed from the event history \eqn{\mathscr{H}_t}; see
#'     \code{\link{redeem_terms}} for available terms.
#'   \item \eqn{\alpha} is the vector of covariate effects.
#'   \item \eqn{\beta_i} and \eqn{\beta_j} are optional actor-specific
#'     baselines (sender and receiver sociality), included via the bare
#'     symbol \code{degree} in the formula.
#'   \item \eqn{f(t, \gamma)} is an optional piecewise-constant temporal
#'     baseline, included via \code{baseline(changepoints)} in the formula.
#' }
#'
#' @details
#' The REM can be viewed as the incidence sub-model of the full
#' \code{\link{dem}}, corresponding to the formation process
#' \eqn{\lambda^{0\rightarrow 1}}. It uses a counting process approach to
#' estimate the influence of various covariates on the timing and occurrence
#' of events, assuming that events are instantaneous points in time.
#'
#' @section Semiparametric Baseline:
#' When \code{semiparametric = TRUE}, the temporal baseline rate of event occurrence
#' is left completely unspecified, and the model parameters are estimated via the Cox
#' partial likelihood using the \code{survival} package. In this path:
#' \itemize{
#'   \item Each observed event time is treated as a failure time, and all non-occurring
#'     dyads at that time constitute the risk set.
#'   \item The exact waiting times between events are conditioned away, meaning that
#'     inference is based solely on the sequence of events and the relative dyadic intensities.
#'   \item This approach is equivalent to the \emph{ordered} (or \emph{conditional})
#'     REM likelihood introduced by Butts (2008). It is highly robust to temporal
#'     fluctuations and baseline misspecification since no piecewise baseline or changepoints
#'     need to be specified.
#'   \item \strong{Limitations}: This path does \strong{not} support the specialized
#'     scalable estimation of sender/receiver popularity effects (\code{degree}) or
#'     piecewise-constant temporal baselines.
#' }
#'
#' @references
#' Fritz, C., Rastelli, R., Fop, M., & Caimo, A. (2026). Scalable Durational Event Models:
#' Application to Physical and Digital Interactions. arXiv:2504.00049.
#'
#' Butts, C. T. (2008). A Relational Event Framework for Social Action.
#' Sociological Methodology, 38(1), 155-200.
#'
#' @param events A matrix of events with columns \code{time}, \code{from},
#'   \code{to}, and optionally \code{type} (1 for start, 3 for exogenous
#'   changes).
#' @param training_start Numeric; the time point at which to start the
#'   estimation. Defaults to 0.
#' @param exogenous_end Numeric; optional end time for exogenous baseline
#'   changes. Defaults to NULL.
#' @param formula A one-sided \code{\link[stats]{formula}} specifying the
#'   sufficient statistics to include in the intensity function. The
#'   right-hand side must be composed of terms from
#'   \code{\link{redeem_terms}}. For example:
#'   \code{~ inertia() + reciprocity() + degree}.
#'   An intercept (\code{~ 1}) is the minimal specification. Defaults to NULL.
#' @param n_nodes Integer; the total number of actors in the network. If \code{NULL} (default), it is automatically identified based on the actors in the \code{events} set.
#' @param directed Logical; whether the interaction events are directed.
#'   Defaults to FALSE.
#' @param coef Numeric vector; initial coefficients for the model.
#'   If provided, this must be a concatenated vector of:
#'   \enumerate{
#'     \item Core coefficients: values for sufficient statistics in the formula.
#'     \item Degree coefficients (if \code{degree} is in the formula): a vector of length \code{n_nodes} (undirected) or \code{2 * n_nodes} (directed, sender effects first then receiver effects).
#'     \item Baseline coefficients (if temporal changepoints are present): a vector of length equal to the number of baseline intervals (equal to number of changepoints if an intercept/degree is present, or changepoints + 1 if neither is present).
#'   }
#'   Defaults to NULL, in which case default starting values are automatically computed.
#' @param semiparametric Logical; whether to use a semiparametric baseline.
#'   Defaults to FALSE. See the 'Semiparametric Baseline' section for details.
#' @param control A list of control parameters from
#'   \code{\link{control.redeem}}. Defaults to \code{control.redeem()}.
#'
#' @return An object of class \code{\link{rem_object}} containing model estimates
#'   and log-likelihoods. See \code{\link{rem_object}} for details on the
#'   components of the returned object and S3 methods.
#'
#' @examples
#' # Simulate some relational event data
#' n <- 20
#' events <- matrix(c(
#'   1.2, 1, 5,
#'   3.1, 2, 8,
#'   4.5, 1, 3
#' ), ncol = 3, byrow = TRUE)
#' colnames(events) <- c("time", "from", "to")
#'
#' # Estimate a simple REM
#' fit <- rem(
#'   events = events,
#'   n_nodes = n,
#'   formula = ~1,
#'   control = control.redeem(it_max = 50)
#' )
#' summary(fit)
#' @importFrom stats formula model.matrix update
#' @export
rem <- function(events,
                training_start = 0,
                exogenous_end = NULL,
                formula = NULL,
                n_nodes = NULL,
                directed = FALSE,
                coef = NULL,
                semiparametric = FALSE,
                control = control.redeem()) {
  call <- match.call()

  processed_actors <- process_event_actors(events, n_nodes, directed = directed)
  events <- processed_actors$events
  n_nodes <- processed_actors$n_nodes

  if (!is.matrix(events)) {
    events <- as.matrix(events)
  }

  # Extract control parameters
  it_max <- control$it_max[1]
  tol <- control$tol[1]
  accelerated <- control$accelerated[1]
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

  if (ncol(events) == 3) {
    events <- cbind(events, 1)
  } else if (ncol(events) != 4) {
    stop("Events matrix must have 3 (time, from, to) or 4 (time, from, to, type) columns.")
  }
  colnames(events) <- c("time", "from", "to", "type")

  if (!is.null(exogenous_end)) {
    if (nrow(events) > 0 && max(events[, 1]) > exogenous_end) {
      warning("exogenous_end is before the last event. Truncating events.")
      events <- events[events[, 1] <= exogenous_end, , drop = FALSE]
    }
    # Add exogenous_end as a pseudo-event early
    events <- rbind(events, matrix(c(exogenous_end, 1, 2, 3), nrow = 1))
    events <- events[order(events[, 1]), ]
  }

  event_numbers <- sum(events[, 4] == 1)

  if (is.null(formula)) {
    return(NULL)
  }

  # Preprocess formulas first to extract changepoints
  preprocessed <- formula_preprocess(
    model_type = "rem",
    formula_1_0 = formula,
    formula_0_1 = formula,
    events = events,
    n_nodes = n_nodes,
    directed = directed
  )

  # Note: in rem(), both formula_0_1 and formula_1_0 in the preprocessor refer to the same formula.
  time_changepoints <- preprocessed$baseline_changepoints_0_1

  # If exogenous changes to the baseline intensity are provided, we add them as artificial events with type 3
  if (length(time_changepoints) > 0) {
    if (estimate == "NR") {
      events <- rbind(events, cbind(time_changepoints, 1, 2, 3))
      events <- events[order(events[, 1]), ]
    }
  }

  time_labels <- preprocessed$baseline_labels_0_1

  if (training_start != 0) {
    if (nrow(events) > 0 && training_start >= max(events[, 1])) {
      stop("The proposed starting time is after the last event.")
    }
    # If we start not from the beginning all change points that are before
    # the training will be deleted
    include <- time_changepoints > training_start
    time_changepoints <- time_changepoints[include]
    if (!is.null(time_labels)) time_labels <- time_labels[include]
  }

  labels_changepoints <- if (!is.null(time_labels)) time_labels else as.character(time_changepoints)
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

  fixed_effects <- (preprocessed$preprocess_0_1$includes_degrees) > 0
  # If the baseline intensity is a constant and no fixed effects are needed,
  # the estimation can be done with glm
  if (!fixed_effects && is.null(time_changepoints)) {
    estimate <- "NR"
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
    simultaneous_interactions = FALSE,
    build_time = build_time,
    max_time = max_time,
    model_type = "rem"
  )

  if (!is.null(build_time) && build_time > 0) {
    preprocessed_tmp <- preprocessed_tmp[preprocessed_tmp$time >= build_time, ]
  }

  if (is.null(preprocessed_tmp) || ncol(preprocessed_tmp) == 0 || nrow(preprocessed_tmp) == 0) {
    stop("Preprocessing failed to produce valid data intervals. Check if your event sequence is valid.")
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
      # More columns than names? Add generic names
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

  if (estimate != "Blockwise" && !is.null(time_changepoints) && length(time_changepoints) > 0) {
    preprocessed_tmp[, time_cat := cut(time_new, breaks = c(-1, time_changepoints, Inf), labels = c("Beg", labels_changepoints))]
    matching <- c(matching, "time_cat")
  }

  preprocessed_tmp <- preprocessed_tmp[!is.na(event), ]

  formula_names <- names(preprocessed$preprocess_0_1$coef_names)
  formula_new <- if (length(formula_names) > 0) {
    formula(paste("~", paste(collapse = "+ ", formula_names)))
  } else {
    ~1
  }
  matching <- c(matching, "pair_id")
  preprocessed_tmp[, weight := 1]
  if (weighting) {
    preprocessed_tmp <- preprocessed_tmp[, .(weight = sum(weight)), by = c(matching)]
  }
  # Prepare data for estimation
  if (estimate != "Blockwise") {
    if (fixed_effects) {
      if (directed) {
        sender_effects <- matrix(0, nrow = nrow(preprocessed_tmp), ncol = n_nodes)
        sender_effects[cbind(seq_len(nrow(sender_effects)), preprocessed_tmp$from)] <- 1
        colnames(sender_effects) <- paste0("sender_", seq_len(n_nodes))

        receiver_effects <- matrix(0, nrow = nrow(preprocessed_tmp), ncol = n_nodes)
        receiver_effects[cbind(seq_len(nrow(receiver_effects)), preprocessed_tmp$to)] <- 1
        colnames(receiver_effects) <- paste0("receiver_", seq_len(n_nodes))

        preprocessed_tmp <- cbind(preprocessed_tmp, sender_effects, receiver_effects)
        formula_new <- update(formula_new, new = paste("~ . -Intercept +", paste0(" sender_", seq_len(n_nodes), " ", collapse = "+"), "+", paste0(" receiver_", seq_len(n_nodes), " ", collapse = "+")))
      } else {
        fixed_effect <- matrix(0, nrow = nrow(preprocessed_tmp), ncol = n_nodes)
        fixed_effect[cbind(rep(seq_len(nrow(fixed_effect)), times = 2), c(preprocessed_tmp$from, preprocessed_tmp$to))] <- 1
        colnames(fixed_effect) <- paste0("effect_", seq_len(n_nodes))
        preprocessed_tmp <- cbind(preprocessed_tmp, fixed_effect)

        formula_new <- update(formula_new, new = paste("~ . -Intercept +", paste0(" effect_", seq_len(n_nodes), " ", collapse = "+")))
      }
    }
    if (!is.null(time_changepoints) && length(time_changepoints) > 0) {
      if (!"time_cat" %in% names(preprocessed_tmp)) {
        preprocessed_tmp$time_cat <- cut(preprocessed_tmp$time_new, breaks = c(-1, time_changepoints, Inf), labels = c("Beg", labels_changepoints))
      }
      preprocessed_tmp <- cbind(preprocessed_tmp, model.matrix(~ -1 + time_cat, preprocessed_tmp))

      # Determine if we need a full baseline (no other non-degree terms)
      full_baseline <- !fixed_effects && length(names(preprocessed$coef_names)) == 0

      if (full_baseline) {
        time_cat_names <- paste0("time_cat", levels(preprocessed_tmp$time_cat))
        formula_new <- formula(paste("~ -1 +", paste(time_cat_names, collapse = " + ")))
      } else {
        time_cat_names <- paste0("time_cat", levels(preprocessed_tmp$time_cat)[-1])
        if (length(time_cat_names) > 0) {
          formula_new <- update(formula_new, new = paste("~ . +", paste0(time_cat_names, collapse = " + ")))
        }
      }
    }
    indicators <- NULL
  } else {
    indicators <- match(attr(terms(formula_new), "term.labels"), names(preprocessed_tmp))
  }

  # Call core estimation helper
  start_time <- Sys.time()
  model <- estimate_transition(
    data = preprocessed_tmp,
    formula_original = formula,
    formula_new = formula_new,
    indicators = indicators,
    n_nodes = n_nodes,
    estimate_method = estimate,
    it_max = it_max,
    tol = tol,
    accelerated = accelerated,
    subsample = subsample,
    verbose = verbose,
    estimate_degree = fixed_effects,
    directed = directed,
    semiparametric = semiparametric,
    labels_changepoints = labels_changepoints,
    time_changepoints = time_changepoints,
    coef_init = coef,
    model_type = "rem",
    process = "0-1",
    return_data = return_data,
    save_hist = save_hist,
    use_glm = use_glm,
    legacy = legacy,
    inf_unidentifiable = inf_unidentifiable,
    events = events
  )
  end_time <- Sys.time()

  res <- list(
    call = call,
    event_numbers = event_numbers,
    model = model,
    events = events,
    formula = formula,
    n_nodes = n_nodes,
    directed = directed,
    build_time = build_time,
    max_time = max_time,
    time_changepoints = time_changepoints,
    labels_changepoints = labels_changepoints,
    training_start = training_start,
    exogenous_end = exogenous_end,
    subsample = subsample,
    return_data = return_data,
    runtime = end_time - start_time,
    window_map = preprocessed$window_map,
    preprocessed = preprocessed
  )

  # Apply descriptive names to results
  if (!is.null(res$model$coefficients)) {
    mapping <- preprocessed$preprocess_0_1$coef_names

    # Use a robust matching for renaming
    curr_names <- names(res$model$coefficients)
    if (!is.null(curr_names)) {
      new_names <- curr_names
      for (i in seq_along(mapping)) {
        internal_name <- names(mapping)[i]
        label <- mapping[i]
        new_names[new_names == internal_name] <- label
      }
      names(res$model$coefficients) <- new_names

      # Cascade renaming to other components
      n_core <- length(res$model$est_core)
      if (n_core > 0) names(res$model$est_core) <- new_names[seq_len(n_core)]

      if (!is.null(res$model$covariance)) {
        n_cov <- nrow(res$model$covariance)
        if (n_cov > 0) {
          colnames(res$model$covariance) <- rownames(res$model$covariance) <- new_names[seq_len(n_cov)]
        }
      }
    }
  }

  class(res) <- "rem"
  res$directed <- directed
  return(res)
}
