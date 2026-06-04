#' @importFrom stats as.formula dpois glm.control glm.fit logLik pnorm terms terms.formula time
#' @importFrom utils tail globalVariables head
#' @importFrom digest digest
#' @importFrom graphics lines
#' @import data.table

# Global variable declarations for data.table
utils::globalVariables(c(
  ".", "weight", "from_to", "time", "event",
  "prediction", "pair_id", "formula_between",
  "formula_tmp", "block_i_nodes", "tmp_type", "head"
))

#' Compute the adjusted rand index (ARI) between two clusterings
#'
#' @noRd
ari <- function(z_star, z) {
  # Compute the contingency table
  cont_table <- table(
    factor(z_star, levels = unique(z_star)),
    factor(z, levels = unique(z))
  )

  # Sum of squares of sums of the contingency table
  sum_comb_c <- sum(choose(colSums(cont_table), 2))
  sum_comb_k <- sum(choose(rowSums(cont_table), 2))
  sum_comb <- sum(choose(cont_table, 2))

  n <- sum(cont_table)
  total_comb <- choose(n, 2)

  # Expected index and max index
  expected_index <- (sum_comb_k * sum_comb_c) / total_comb
  max_index <- (sum_comb_k + sum_comb_c) / 2

  # Adjusted Rand Index
  ari <- (sum_comb - expected_index) / (max_index - expected_index)
  return(ari)
}


.split_plus <- function(expr) {
  out <- list()
  rec <- function(e) {
    if (is.call(e) && identical(e[[1L]], as.name("+"))) {
      rec(e[[2L]])
      rec(e[[3L]])
    } else {
      out[[length(out) + 1L]] <<- e
    }
  }
  rec(expr)
  out
}

.deparse1 <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = "")

#' Preprocess a single formula for model terms
#'
#' @description
#' This function takes an R formula and extracts the necessary information to build the model's
#' design matrix. It identifies special terms, transformations, and associated data.
#'
#' @param formula An R formula object.
#' @param n_nodes Number of nodes.
#' @param model_type Either "dem" or "rem".
#' @param process Either "0-1" (incidence) or "1-0" (duration).
#' @param directed Logical; if TRUE, the model is directed (defaults to FALSE).
#' @param simulation Logical; if TRUE, the formula is being preprocessed for simulation (defaults to FALSE).
#' @keywords internal
formula_preprocess_single <- function(formula, n_nodes, model_type = "dem", process = "0-1", directed = FALSE, simulation = FALSE) {
  time_changepoints <- NULL
  formula <- stats::as.formula(formula)
  term_labels <- attr(stats::terms(formula), "term.labels")
  has_degree <- "degree" %in% term_labels
  has_degrees <- "degrees" %in% term_labels

  if (has_degree || has_degrees) {
    if (has_degree) formula <- stats::update(formula, . ~ . - degree)
    if (has_degrees) formula <- stats::update(formula, . ~ . - degrees)
    includes_degrees <- TRUE
  } else {
    includes_degrees <- FALSE
  }
  env_eval <- environment(formula)
  if (is.null(env_eval)) env_eval <- parent.frame()

  formula_info <- rhs_terms_as_list(formula, n_nodes = n_nodes, model_type = model_type, process = process, directed = directed, env = env_eval)

  # Detect changepoints from covariates
  cov_changepoints <- NULL
  for (i in seq_along(formula_info)) {
    if (is.list(formula_info[[i]]$data)) {
      nms <- names(formula_info[[i]]$data)
      if (!is.null(nms)) {
        suppressWarnings(ts <- as.numeric(nms))
        ts <- ts[!is.na(ts) & ts > 1e-10]
        cov_changepoints <- c(cov_changepoints, ts)
      }
    }
  }

  # Extract baseline changepoints and labels if present (robust check for base_name and presence of changepoints)
  is_baseline_term <- sapply(formula_info, function(x) identical(x$base_name, "baseline") && !is.null(x[["changepoints"]]))
  has_baseline_term <- any(is_baseline_term) || any(sapply(formula_info, function(x) identical(x$base_name, "baseline")))

  time_labels <- NULL
  if (any(is_baseline_term)) {
    baseline_info <- formula_info[[which(is_baseline_term)[1]]]
    time_changepoints <- baseline_info$changepoints
    time_labels <- baseline_info$labels
    formula_info <- formula_info[!is_baseline_term]
  }

  has_any_changepoints <- !is.null(time_changepoints) || length(cov_changepoints) > 0

  # NEVER allow Intercept and fixed effects in the model (except for simulation)
  has_explicit_intercept <- any(sapply(formula_info, function(x) tolower(x$base_name) == "intercept"))
  wants_intercept <- attr(stats::terms(formula), "intercept") == 1
  intercept_removed <- FALSE
  if ((includes_degrees || has_baseline_term) && !simulation) {
    if (has_explicit_intercept) {
      formula_info <- formula_info[!sapply(formula_info, function(x) tolower(x$base_name) == "intercept")]
      has_explicit_intercept <- FALSE
      intercept_removed <- TRUE
    } else if (wants_intercept) {
      intercept_removed <- TRUE
    }
  }

  # Prepend Intercept if it's missing (only for simulation)
  append_intercept <- FALSE
  if (!has_explicit_intercept && wants_intercept && simulation && (has_baseline_term || has_any_changepoints)) {
    append_intercept <- TRUE
  }

  if (append_intercept) {
    arg_vals <- list(label = "Intercept", base_name = "Intercept")
    attr(arg_vals, "process") <- process
    intercept_term <- InitRedeemTerm("Intercept", arg_vals,
      model_type = model_type, process = process, n_nodes = n_nodes, directed = directed
    )
    intercept_term$label <- "Intercept"
    formula_info <- c(list(Intercept = intercept_term), formula_info)
    has_explicit_intercept <- TRUE
  }

  # ALWAYS move Intercept to the front if present
  if (has_explicit_intercept) {
    int_idx <- which(sapply(formula_info, function(x) tolower(x$base_name) == "intercept"))
    if (length(int_idx) > 0 && int_idx[1] != 1) {
      formula_info <- c(formula_info[int_idx[1]], formula_info[-int_idx[1]])
    }
  }


  transformation_list <- lapply(formula_info, function(x) x$transformation)
  transformation_list <- unlist(lapply(transformation_list, function(x) {
    if (is.null(x)) {
      "identity"
    } else {
      x
    }
  }))

  term_per_term <- sapply(formula_info, function(x) {
    if (!is.null(x$window) && is.finite(x$window)) paste0(x$base_name, "_wt", x$window) else x$base_name
  })
  type_per_term <- lapply(formula_info, function(x) x$type)
  type_per_term <- unlist(lapply(type_per_term, function(x) {
    if (is.null(x)) {
      1
    } else {
      x
    }
  }))
  window_per_term <- unlist(lapply(formula_info, function(x) {
    if (is.null(x$window)) Inf else x$window
  }))

  data_per_term <- lapply(formula_info, function(x) {
    res <- if (is.list(x$data)) x$data else if (!is.null(x$eval_at_zero)) x$eval_at_zero else if (!is.null(x$data)) x$data else matrix(1)
    if (is.list(res)) {
      # For list data (time-varying), we keep it as a list but ensure each element is a matrix
      res <- lapply(res, function(m) {
        if (!is.matrix(m)) m <- as.matrix(m)
        if (!is.numeric(m)) storage.mode(m) <- "double"
        return(m)
      })
    } else {
      if (!is.matrix(res)) res <- as.matrix(res)
      if (!is.numeric(res)) storage.mode(res) <- "double"
    }
    return(res)
  })
  name_per_term <- as.character(names(formula_info))
  # names(name_per_term) are the internal names (e.g. dyadic_cov_identity)
  # we ensure they are slightly sanitized for data.frame/matching use,
  # but they are NOT the ones shown in summary.
  if (!is.null(name_per_term)) {
    name_per_term <- gsub(pattern = "(", replacement = ".", x = name_per_term, fixed = TRUE)
    name_per_term <- gsub(pattern = ")", replacement = ".", x = name_per_term, fixed = TRUE)
    name_per_term <- gsub(pattern = "=", replacement = ".", x = name_per_term, fixed = TRUE)
    name_per_term <- gsub(pattern = " ", replacement = ".", x = name_per_term, fixed = TRUE)
    # Ensure values and names are the same
    names(name_per_term) <- name_per_term
  }
  term_per_term <- as.character(term_per_term)

  stream_per_term <- lapply(formula_info, function(x) x$event_stream)

  return(list(
    data_list = data_per_term,
    transformation_list = as.character(unlist(transformation_list)),
    type_list = type_per_term,
    coef_names = name_per_term,
    term_names = term_per_term,
    window_list = window_per_term,
    stream_per_term = stream_per_term,
    includes_degrees = includes_degrees,
    intercept_removed = intercept_removed,
    time_changepoints = time_changepoints,
    time_labels = time_labels,
    has_baseline_term = has_baseline_term
  ))
}


#' Convert right-hand side of a formula to a list of term information
#'
#' @param formula The formula to parse.
#' @param env The environment in which to evaluate terms.
#' @param evaluate_calls Logical; if `TRUE`, evaluates the full calls.
#' @param model_type Either "dem" or "rem".
#' @param process Either "0-1" or "1-0".
#' @param directed Logical; if TRUE, the model is directed.
#'
#' @param n_nodes Number of nodes.
#'
#' @return A list of term information, including labels, base names, and data.
#' @keywords internal
rhs_terms_as_list <- function(formula, n_nodes, env = NULL, evaluate_calls = FALSE, model_type = "dem", process = "0-1", directed = FALSE) {
  formula <- as.formula(formula)
  if (is.null(env)) {
    env <- environment(formula)
    if (is.null(env)) env <- parent.frame()
  }
  rhs_expr <- if (length(formula) == 3) formula[[3L]] else formula[[2L]]
  terms_exprs <- .split_plus(rhs_expr)

  out <- list()
  taken_names <- character(0L)

  for (term_expr in terms_exprs) {
    if (is.numeric(term_expr) && term_expr == 1) {
      arg_vals <- list(label = "Intercept", base_name = "Intercept")
      attr(arg_vals, "process") <- process
      out[["Intercept"]] <- InitRedeemTerm("Intercept", arg_vals,
        model_type = model_type, process = process, n_nodes = n_nodes, directed = directed
      )
      out[["Intercept"]]$label <- "Intercept"
      next
    }
    if (is.symbol(term_expr)) {
      # cat("Symbol")
      # Skip purely numerical constants (like -1 or 0)
      if (is.numeric(term_expr)) {
        next
      }
      if (is.call(term_expr) && length(term_expr) == 2 && is.numeric(term_expr[[2L]])) {
        next
      }
      base_name <- as.character(term_expr)
      if (tolower(base_name) == "intercept") {
        arg_vals <- list(label = "Intercept", base_name = "Intercept")
        attr(arg_vals, "process") <- process
        out[["Intercept"]] <- InitRedeemTerm("Intercept", arg_vals,
          model_type = model_type, process = process, n_nodes = n_nodes, directed = directed
        )
        out[["Intercept"]]$label <- "Intercept"
        next
      }
      if (base_name %in% taken_names) {
        next
      }
      taken_names <- c(taken_names, base_name)

      v <- try(eval(term_expr, envir = env), silent = TRUE)

      # For symbols, we still want to go through the InitRedeemTerm logic if possible
      arg_vals <- list(base_name = base_name, label = .deparse1(term_expr))
      if (!inherits(v, "try-error") && !is.null(v)) {
        arg_vals$data <- v
      }
      attr(arg_vals, "process") <- process

      init_res <- try(InitRedeemTerm(base_name, arg_vals, model_type = model_type, process = process, n_nodes = n_nodes, directed = directed), silent = TRUE)

      if (inherits(init_res, "try-error")) {
        # If no Init function, use default behavior
        out[[base_name]] <- list(
          label = .deparse1(term_expr),
          base_name = base_name,
          data = if (inherits(v, "try-error") || is.null(v)) matrix(1) else v,
          eval_at_zero = if (inherits(v, "try-error") || is.null(v)) matrix(1, n_nodes, n_nodes) else v
        )
      } else {
        out[[base_name]] <- c(list(label = .deparse1(term_expr)), init_res)
      }
    } else if (is.call(term_expr)) {
      # cat("call term")
      fun_sym <- term_expr[[1L]]
      base_name <- if (is.symbol(fun_sym)) as.character(fun_sym) else .deparse1(fun_sym)

      if (tolower(base_name) == "intercept") {
        arg_vals <- list(label = "Intercept", base_name = "Intercept")
        attr(arg_vals, "process") <- process
        out[["Intercept"]] <- InitRedeemTerm("Intercept", arg_vals,
          model_type = model_type, process = process, n_nodes = n_nodes, directed = directed
        )
        out[["Intercept"]]$label <- "Intercept"
        next
      }

      # raw argument expressions
      arg_exprs <- as.list(term_expr)[-1L]
      arg_names <- names(arg_exprs)
      if (is.null(arg_names)) arg_names <- rep("", length(arg_exprs))

      # Detect if the 'data' argument is a symbol to use in the term name
      data_arg_symbol <- NULL
      data_idx <- which(arg_names == "data")
      if (length(data_idx) == 0 && length(arg_exprs) > 0 && arg_names[1] == "") {
        data_idx <- 1
      }
      if (length(data_idx) > 0 && is.symbol(arg_exprs[[data_idx]])) {
        data_arg_symbol <- as.character(arg_exprs[[data_idx]])
      }

      # prepare container with named, evaluated arguments
      # unnamed arguments get positional names data, ..2, ...
      pos_names <- ifelse(arg_names == "", ifelse(seq_along(arg_exprs) == 1, "data", paste0("..", seq_along(arg_exprs))), arg_names)
      arg_vals <- vector("list", length(arg_exprs))
      names(arg_vals) <- pos_names

      for (i in seq_along(arg_exprs)) {
        v <- try(eval(arg_exprs[[i]], envir = env), silent = TRUE)
        arg_vals[[i]] <- if (inherits(v, "try-error")) NULL else v
      }

      # Add metadata for InitRedeemTerm
      arg_vals$base_name <- base_name
      arg_vals$label <- gsub(pattern = '\\\"', replacement = "'", x = .deparse1(term_expr))
      attr(arg_vals, "process") <- process

      init_res <- InitRedeemTerm(base_name, arg_vals, model_type = model_type, process = process, n_nodes = n_nodes, directed = directed)

      # optionally evaluate the whole call (if it's not a redeem_term function)
      evaluated <- NULL
      if (evaluate_calls) {
        tmp <- try(eval(term_expr, envir = env), silent = TRUE)
        if (!inherits(tmp, "try-error")) evaluated <- tmp
      }

      # Refine base_name if data argument is a symbol
      actual_base_name <- if (!is.null(init_res$base_name)) init_res$base_name else base_name
      if (!is.null(data_arg_symbol) && actual_base_name %in% c("dyadic_cov", "monadic_cov", "node_cov", "interaction_cov")) {
        actual_base_name <- paste0(actual_base_name, "_", data_arg_symbol)
        init_res$base_name <- actual_base_name
      }


      entry <- c(
        label = arg_vals$label,
        init_res
      )
      if (!is.null(evaluated)) entry$.evaluated <- evaluated

      # Determine unique name for the list entry
      name_addon <- ""
      if (!is.null(init_res$type)) {
        name_addon <- paste0(name_addon, init_res$type)
      }
      if (!is.null(init_res$window) && is.finite(init_res$window)) {
        name_addon <- paste0(name_addon, ifelse(name_addon == "", "", "_"), "wt", init_res$window)
      }
      if (!is.null(init_res$data) && is.character(init_res$data)) {
        name_addon <- paste0(name_addon, ifelse(name_addon == "", "", "_"), init_res$data)
      }
      if (!is.null(init_res$transformation)) {
        name_addon <- paste0(name_addon, ifelse(name_addon == "", "", "_"), init_res$transformation)
      }
      if (!is.null(init_res$variant)) {
        name_addon <- paste0(name_addon, ifelse(name_addon == "", "", "_"), init_res$variant)
      }

      actual_base_name <- if (!is.null(init_res$base_name)) init_res$base_name else base_name
      if (name_addon != "") {
        base_elt_name <- paste0(actual_base_name, "_", name_addon)
      } else {
        # Fallback to label hash if no other identifiers are available to avoid collisions
        if (actual_base_name %in% c("dyadic_cov", "node_cov", "interaction_cov") ||
          grepl("^(dyadic_cov|node_cov|interaction_cov|monadic_cov)_", actual_base_name)) {
          hash <- substr(digest::digest(arg_vals$label), 1, 6)
          base_elt_name <- paste0(actual_base_name, "_", hash)
        } else {
          base_elt_name <- actual_base_name
        }
      }
      elt_name <- base_elt_name
      count <- 1
      while (elt_name %in% taken_names) {
        elt_name <- paste0(base_elt_name, ".", count)
        count <- count + 1
      }
      taken_names <- c(taken_names, elt_name)
      out[[elt_name]] <- entry
    }
  }
  class(out) <- "redeem.formulainfo"
  out
}

#' Preprocess Formulas for Model Terms with Event and Node Information
#'
#' This function processes two model formulas, each of which can specify transformation
#' and data arguments, and combines the preprocessed results with additional event and node information.
#' It identifies unique coefficient names across both formulas, determining which terms to include
#' based on uniqueness, and returns structured lists for data, transformations, and term names.
#'
#' @param formula_0_1 An R formula for the `0 -> 1` model terms, also specifying terms with
#'        optional `transformation = ...` and `data = ...` arguments.
#' @param formula_1_0 An R formula specifying model terms, with the form `~ term1 + term2 + ...`,
#'        including optional `transformation = ...` and `data = ...` arguments. Typically, this formula
#'        represents the model's `1 -> 0` terms.
#' @param events A data frame or list representing the events associated with the model,
#'        for use in estimating or evaluating the model.
#' @param n_nodes An integer specifying the number of nodes in the network or model structure.
#' @param exo_breaks Optional; a vector or list specifying external breaks, if applicable,
#'        for segmenting the data across the terms.
#'
#' @return A list containing the following components:
#' \item{events}{The events input, retained for use in model estimation or evaluation.}
#' \item{n_nodes}{The number of nodes specified in the input.}
#' \item{data_list}{A combined list of matrices for each term’s data, from both formulas,
#' where each matrix corresponds to the data for a specific term.}
#' \item{transformation_list}{A combined character vector of transformation types for
#' each term, with `"identity"` for terms without specified transformations.}
#' \item{coef_names}{A character vector of coefficient names for each term, combining terms
#' across both formulas and ensuring uniqueness.}
#' \item{term_names}{A character vector of term names, ordered to match `data_list`.}
#' \item{preprocess_1_0}{The output list from `formula_preprocess_single` applied to `formula_1_0`.}
#' \item{preprocess_0_1}{The output list from `formula_preprocess_single` applied to `formula_0_1`.}
#' \item{included_1_0}{A logical vector indicating whether each term in `coef_names`
#' comes from `formula_1_0`.}
#' \item{included_0_1}{A logical vector indicating whether each term in `coef_names`
#' comes from `formula_0_1`.}
#'
#' @details
#' The function first calls `formula_preprocess_single` on `formula_1_0` and `formula_0_1`
#' separately to obtain individual term processing details. It then identifies unique terms
#' across both formulas and combines the term data, transformations, and coefficient names
#' into a single output list, structured for use in further modeling or evaluation.
#'
#' @examples
#' \dontrun{
#' # Define simple event data
#' event_data <- matrix(c(
#'   1.2, 1, 5, 1,
#'   2.5, 1, 5, 0,
#'   3.1, 2, 8, 1,
#'   4.4, 2, 8, 0
#' ), ncol = 4, byrow = TRUE)
#' colnames(event_data) <- c("time", "from", "to", "type")
#'
#' # Preprocess the formulas
#' formula_preprocess(
#'   formula_1_0 = ~ current_interaction() + current_common_partners(),
#'   formula_0_1 = ~ general_common_partners(),
#'   events = event_data,
#'   n_nodes = 10
#' )
#' }
#'
#' @seealso \code{\link{formula_preprocess_single}}
#' @param formula_0_1 Optional; an R formula for the `0 -> 1` terms.
#' @param model_type Either "dem" or "rem".
#' @param formula_1_0 Optional; an R formula for the `1 -> 0` terms.
#' @param events A data frame or list representing the events.
#' @param n_nodes An integer specifying the number of nodes.
#' @param exo_breaks Optional; a vector or list specifying external breaks.
#' @param directed Logical; if TRUE, the model is directed (defaults to FALSE).
#' @param simulation Logical; if TRUE, the formula is being preprocessed for simulation (defaults to FALSE).
#' @keywords internal
#' @return A list containing the preprocessed information.
formula_preprocess <- function(formula_0_1 = NULL,
                               model_type = "dem",
                               formula_1_0 = NULL,
                               events = matrix(c(0, 0, 0), nrow = 1),
                               n_nodes,
                               exo_breaks = NULL,
                               directed = FALSE,
                               simulation = FALSE) {
  time_changepoints <- NULL
  if (is.null(formula_0_1) && is.null(formula_1_0)) {
    stop("At least one formula (formula_0_1 or formula_1_0) must be provided.")
  }

  preprocess_1_0 <- if (!is.null(formula_1_0)) {
    formula_preprocess_single(formula_1_0, n_nodes = n_nodes, model_type = model_type, process = "1-0", directed = directed, simulation = simulation)
  } else {
    NULL
  }

  preprocess_0_1 <- if (!is.null(formula_0_1)) {
    formula_preprocess_single(formula_0_1, n_nodes = n_nodes, model_type = model_type, process = "0-1", directed = directed, simulation = simulation)
  } else {
    NULL
  }


  fixed_effects <- (isTRUE(preprocess_0_1$includes_degrees) ||
    isTRUE(preprocess_1_0$includes_degrees))

  has_baseline <- isTRUE(preprocess_0_1$has_baseline_term) || isTRUE(preprocess_1_0$has_baseline_term)

  if ((fixed_effects || has_baseline) && !simulation) {
    if (!is.null(formula_1_0)) formula_1_0 <- update(formula_1_0, ~ . - Intercept)
    if (!is.null(formula_0_1)) formula_0_1 <- update(formula_0_1, ~ . - Intercept)
  }

  all_included <- c(preprocess_1_0$coef_names, preprocess_0_1$coef_names)
  all_unique <- unique(all_included)
  if ("Intercept" %in% all_unique) {
    all_unique <- c("Intercept", setdiff(all_unique, "Intercept"))
  }
  if (is.null(all_unique)) all_unique <- character(0)

  if (!is.null(preprocess_1_0)) {
    preprocess_1_0$include <- preprocess_1_0$coef_names %in% all_unique
    all_unique <- all_unique[!all_unique %in% preprocess_1_0$coef_names]
  }

  if (!is.null(preprocess_0_1)) {
    preprocess_0_1$include <- preprocess_0_1$coef_names %in% all_unique
  }

  res <- list(
    events = events,
    n_nodes = n_nodes,
    data_list = c(preprocess_1_0$data_list[preprocess_1_0$include], preprocess_0_1$data_list[preprocess_0_1$include]),
    transformation_list = c(
      preprocess_1_0$transformation_list[preprocess_1_0$include],
      preprocess_0_1$transformation_list[preprocess_0_1$include]
    ),
    coef_names = c(preprocess_1_0$coef_names[preprocess_1_0$include], preprocess_0_1$coef_names[preprocess_0_1$include]),
    term_names = c(preprocess_1_0$term_names[preprocess_1_0$include], preprocess_0_1$term_names[preprocess_0_1$include]),
    stream_list = c(preprocess_1_0$stream_per_term[preprocess_1_0$include], preprocess_0_1$stream_per_term[preprocess_0_1$include]),
    preprocess_1_0 = preprocess_1_0,
    preprocess_0_1 = preprocess_0_1,
    baseline_changepoints_0_1 = preprocess_0_1$time_changepoints,
    baseline_changepoints_1_0 = preprocess_1_0$time_changepoints,
    baseline_labels_0_1 = preprocess_0_1$time_labels,
    baseline_labels_1_0 = preprocess_1_0$time_labels
  )

  # Windowed statistics: inject dissolution events into all relevant streams
  all_windows <- c(
    preprocess_1_0$window_list[preprocess_1_0$include],
    preprocess_0_1$window_list[preprocess_0_1$include]
  )
  unique_finite_windows <- unique(all_windows[is.finite(all_windows)])

  if (length(unique_finite_windows) > 0) {
    # Assign unique type to each window (starting from 10)
    window_map <- stats::setNames(as.numeric(10 + seq_along(unique_finite_windows)), unique_finite_windows)
    res$window_map <- window_map

    # Function to inject dissolution events into a stream
    global_max_time <- if (nrow(events) > 0) max(events[, 1], na.rm = TRUE) else 0
    inject_dissolutions <- function(stream) {
      if (is.null(stream) || nrow(stream) == 0) {
        return(stream)
      }

      # Use column 4 as event type (1=formation)
      if (ncol(stream) < 4) {
        return(stream)
      }
      formation_events <- stream[stream[, 4] == 1, , drop = FALSE]
      if (nrow(formation_events) == 0) {
        return(stream)
      }

      diss_list <- list()
      for (w in unique_finite_windows) {
        w_type <- window_map[as.character(w)]
        new_events <- formation_events
        new_events[, 1] <- new_events[, 1] + w
        # Add +10 to w_type to distinguish windowed expirations from formation (1) and dissolution (0)
        new_events[, 4] <- w_type + 10

        # Constraint: must be <= global_max_time
        new_events <- new_events[new_events[, 1] <= global_max_time, , drop = FALSE]
        if (nrow(new_events) > 0) diss_list[[as.character(w)]] <- new_events
      }

      if (length(diss_list) > 0) {
        stream <- rbind(stream, do.call(rbind, diss_list))
        stream <- stream[order(stream[, 1], stream[, 4]), , drop = FALSE]
      }
      return(stream)
    }

    # Inject into main events
    events <- inject_dissolutions(events)

    # Inject into term-specific streams
    if (!is.null(res$stream_list)) {
      res$stream_list <- lapply(res$stream_list, inject_dissolutions)
    }
  } else {
    res$window_map <- numeric(0)
  }

  res$term_names <- c(
    res$preprocess_1_0$term_names[res$preprocess_1_0$include],
    res$preprocess_0_1$term_names[res$preprocess_0_1$include]
  )

  # Ensure character vectors are not NULL for Rcpp
  if (is.null(res$coef_names)) res$coef_names <- character(0)
  if (is.null(res$term_names)) res$term_names <- character(0)
  if (is.null(res$transformation_list)) res$transformation_list <- character(0)
  if (is.null(res$stream_list)) res$stream_list <- vector("list", 0)

  if (length(res$data_list) > 0) {
    if (FALSE %in% unique(unlist(lapply(res$data_list, dim))) %in% (c(1, n_nodes))) {
      stop("Some of the covariate data is of the wrong format! Please check.")
    }
  }

  res$events <- events
  res$included_1_0 <- if (!is.null(res$coef_names)) res$coef_names %in% res$preprocess_1_0$coef_names else logical(0)
  res$included_0_1 <- if (!is.null(res$coef_names)) res$coef_names %in% res$preprocess_0_1$coef_names else logical(0)
  res$formula_0_1 <- formula_0_1
  res$formula_1_0 <- formula_1_0
  return(res)
}

get_start_coefs <- function(formula,
                            n_nodes,
                            fixed_effects,
                            time_changepoints,
                            directed = FALSE,
                            model_type = "dem",
                            process = "0-1",
                            simulation = FALSE) {
  term_labels <- attr(stats::terms(formula), "term.labels")
  if ("degree" %in% term_labels) formula <- stats::update(formula, . ~ . - degree)
  if ("degrees" %in% term_labels) formula <- stats::update(formula, . ~ . - degrees)

  # Use the same logic as formula_preprocess_single to identify terms
  formula_info <- rhs_terms_as_list(formula, n_nodes = n_nodes, model_type = model_type, process = process, directed = directed)

  # Extract baseline info (robust check for presence of changepoints)
  is_baseline_term <- sapply(formula_info, function(x) identical(x$base_name, "baseline") && !is.null(x[["changepoints"]]))
  has_baseline_term <- any(is_baseline_term) || any(sapply(formula_info, function(x) identical(x$base_name, "baseline")))
  n_baseline <- sum(is_baseline_term)
  formula_info <- formula_info[!is_baseline_term]

  # Detect changepoints from covariates
  cov_changepoints <- NULL
  for (i in seq_along(formula_info)) {
    if (is.list(formula_info[[i]]$data)) {
      nms <- names(formula_info[[i]]$data)
      if (!is.null(nms)) {
        suppressWarnings(ts <- as.numeric(nms))
        ts <- ts[!is.na(ts) & ts > 1e-10]
        cov_changepoints <- c(cov_changepoints, ts)
      }
    }
  }
  has_any_changepoints <- !is.null(time_changepoints) || length(cov_changepoints) > 0

  # Kick out Intercept if fixed effects or baseline are present (matching formula_preprocess_single logic)
  has_explicit_intercept <- any(sapply(formula_info, function(x) tolower(x$base_name) == "intercept"))
  if ((fixed_effects || has_baseline_term) && has_explicit_intercept && !simulation) {
    formula_info <- formula_info[!sapply(formula_info, function(x) tolower(x$base_name) == "intercept")]
    has_explicit_intercept <- FALSE
  }

  # Prepend Intercept if it's missing (only for simulation)
  wants_intercept <- attr(stats::terms(formula), "intercept") == 1
  append_intercept <- FALSE
  if (!has_explicit_intercept && wants_intercept && simulation && (has_baseline_term || has_any_changepoints)) {
    append_intercept <- TRUE
  }
  if (append_intercept) {
    n_coef_core <- length(formula_info) + 1
  } else {
    n_coef_core <- length(formula_info)
  }

  coef_core <- rep(0, n_coef_core)

  if (fixed_effects) {
    coef_degree <- rep(-5, if (directed) 2 * n_nodes else n_nodes)
  } else {
    coef_degree <- numeric(0)
  }

  if (has_any_changepoints) {
    # If there's no intercept/fixed effects, it's a full baseline.
    full_baseline <- !(fixed_effects || (has_explicit_intercept && !simulation) || append_intercept)
    n_time <- length(time_changepoints) + (if (full_baseline) 1 else 0)
    coef_time <- rep(0, n_time)
  } else {
    coef_time <- numeric(0)
  }

  return(list(
    coef_core = coef_core,
    coef_degree = coef_degree,
    coef_time = coef_time
  ))
}


#' Validate the Structure of a Durational Event List
#'
#' This function checks the validity of a dyadic interaction matrix by ensuring
#' that each interaction start event has a corresponding end event, that no
#' interactions overlap within a dyad, and that no missing values are present.
#'
#' @param df A data frame with at least four columns, representing events, where:
#' \describe{
#'   \item{Column 1}{Event timing or ID (e.g., timestamp).}
#'   \item{Column 2}{"From" node ID for the dyadic interaction.}
#'   \item{Column 3}{"To" node ID for the dyadic interaction.}
#'   \item{Column 4}{Event type (1 for start, 0 for end of interaction).}
#' }
#'
#' @return Logical; \code{TRUE} if all interactions are valid, \code{FALSE} otherwise.
#' If the data contains missing values, the function issues a warning and returns \code{FALSE}.
#'
#' @details
#' The function performs the following checks:
#' \itemize{
#'   \item Missing values: If any are found, a warning is issued and \code{FALSE} is returned.
#'   \item Interaction pairing: Each start event (1) must have a corresponding end event (0) without overlap.
#'   \item Non-overlapping intervals: Ensures that no start event occurs while another interaction is active.
#' }
#'
#' @param return_matrix Logical; if TRUE, returns the (potentially repaired) event matrix.
#'   Defaults to FALSE, in which case the function returns \code{TRUE} if the matrix is valid.
#' @param start_time Numeric; optional reference time for adding missing start events.
#'   Defaults to NULL, in which case the earliest time in the data is used.
#' @examples
#' # Create a valid event matrix with durational events (start=1, end=0)
#' df <- matrix(c(
#'   1.0, 1, 2, 1,
#'   2.0, 1, 2, 0,
#'   1.5, 3, 4, 1,
#'   3.0, 3, 4, 0
#' ), ncol = 4, byrow = TRUE)
#' colnames(df) <- c("time", "from", "to", "type")
#'
#' # Check if the event matrix is valid
#' check_matrix(df)
#' @export
check_matrix <- function(df, return_matrix = FALSE, start_time = NULL) {
  # Convert to matrix if it's a data.frame or data.table to ensure consistent indexing
  if (!is.matrix(df)) df <- as.matrix(df)

  if (ncol(df) < 4) {
    stop(sprintf("Event matrix must have at least 4 columns (time, from, to, type). Found only %d columns.", ncol(df)))
  }

  # Detect columns by name if possible, otherwise assume standard order (1=time, 2=from, 3=to, 4=type)
  col_names <- tolower(colnames(df))
  time_col <- if ("time" %in% col_names) which(col_names == "time")[1] else 1
  from_col <- if ("from" %in% col_names) which(col_names == "from")[1] else 2
  to_col <- if ("to" %in% col_names) which(col_names == "to")[1] else 3
  type_col <- if (any(col_names %in% c("type", "event"))) {
    which(col_names %in% c("type", "event"))[1]
  } else {
    4
  }

  # Sort by From, To, then Time, then Type (1 before 0) to ensure stable processing
  # We use as.numeric() for columns to handle cases where df is a character matrix.
  df <- df[order(df[, from_col], df[, to_col], as.numeric(df[, time_col]), -as.numeric(df[, type_col])), , drop = FALSE]

  errors <- character()

  if (any(is.na(df))) {
    na_rows <- which(apply(df, 1, function(x) any(is.na(x))))
    errors <- c(errors, sprintf(
      "Found missing values in the data at row(s): %s",
      paste(head(na_rows, 10), collapse = ", ")
    ))
  }

  if (is.null(start_time)) {
    # Ensure we use the detected time column
    start_time <- if (nrow(df) > 0) min(as.numeric(df[, time_col]), na.rm = TRUE) else 0
  }

  # Get the unique pairs of 'from' and 'to'
  unique_dyads <- unique(df[, c(from_col, to_col), drop = FALSE])
  added_starts <- list()
  warnings_list <- character()

  # Track indices to keep
  keep_idx <- rep(TRUE, nrow(df))

  for (dyad_idx in seq_len(nrow(unique_dyads))) {
    from <- unique_dyads[dyad_idx, 1]
    to <- unique_dyads[dyad_idx, 2]

    # Indices for this dyad
    dyad_indices <- which(df[, from_col] == from & df[, to_col] == to)
    dyad_data <- df[dyad_indices, , drop = FALSE]

    # Subset to interaction events (0 or 1)
    interaction_mask <- as.numeric(dyad_data[, type_col]) %in% c(0, 1)
    interaction_indices <- dyad_indices[interaction_mask]
    interaction_events <- dyad_data[interaction_mask, , drop = FALSE]

    open_interaction <- FALSE # Track if there's an ongoing interaction
    last_start_idx <- NA
    last_start_time <- NA

    for (i in seq_len(nrow(interaction_events))) {
      time <- as.numeric(interaction_events[i, time_col])
      event_type <- as.numeric(interaction_events[i, type_col])
      current_idx <- interaction_indices[i]

      if (event_type == 1) {
        # Start of an interaction
        if (open_interaction) {
          # If there's already an open interaction, this "start" would overlap
          errors <- c(errors, sprintf(
            "Overlap found for dyad (%s, %s): Start event at t=%.4f while interaction from t=%.4f is still open.",
            as.character(from), as.character(to), as.numeric(time), as.numeric(last_start_time)
          ))
        }
        open_interaction <- TRUE
        last_start_idx <- current_idx
        last_start_time <- time
      } else if (event_type == 0) {
        # End of an interaction
        if (!open_interaction) {
          # REPAIR (SILENTLY): Add a start event at the observation start or the end event time
          # No warning is issued for these cases as requested.
          repair_time <- min(start_time, time)
          added_starts[[length(added_starts) + 1]] <- c(repair_time, from, to, 1)
          open_interaction <- TRUE
        }
        # Close the open interaction
        open_interaction <- FALSE
      }
    }

    # Unclosed interactions are kept as censored data silently (as requested)
    if (open_interaction) {
      # No action needed - just leave open_interaction as TRUE
    }
  }

  if (length(added_starts) > 0) {
    # Combine with original data
    new_starts_mat <- matrix(NA, nrow = length(added_starts), ncol = ncol(df))
    for (i in seq_along(added_starts)) {
      new_starts_mat[i, time_col] <- added_starts[[i]][1]
      new_starts_mat[i, from_col] <- added_starts[[i]][2]
      new_starts_mat[i, to_col] <- added_starts[[i]][3]
      new_starts_mat[i, type_col] <- added_starts[[i]][4]
    }
    # Ensure column names match if any
    colnames(new_starts_mat) <- colnames(df)

    # Final idempotency filter: don't add rows that already exist
    df <- rbind(df, new_starts_mat)
    df <- df[!duplicated(df[, c(time_col, from_col, to_col, type_col)]), , drop = FALSE]
  }

  # Final sort using Time as the primary key, then Type (1 before 0) for stability
  # This ensures the matrix is chronologically ordered for dem() and rem().
  df <- df[order(as.numeric(df[, time_col]), -as.numeric(df[, type_col]), df[, from_col], df[, to_col]), , drop = FALSE]

  # Issue errors only (warnings for repairs/unclosed are now suppressed)
  for (err in unique(errors)) warning(err)

  if (return_matrix) {
    if (length(errors) > 0) attr(df, "has_errors") <- TRUE
    return(df)
  }

  if (length(errors) > 0) {
    return(FALSE)
  }

  return(TRUE)
}

# Union Dyadic Datasets
union_dyadic_datasets <- function(dt1, dt2,
                                  topo_keys = c("pair_id", "from", "to", "time", "time_new", "event", "status", "from_avail", "to_avail"),
                                  suffix1 = "d1", suffix2 = "d2",
                                  impute_zero = FALSE, omit_na = FALSE) {
  # 1. Validation & Safety
  if (!data.table::is.data.table(dt1)) dt1 <- data.table::as.data.table(dt1)
  if (!data.table::is.data.table(dt2)) dt2 <- data.table::as.data.table(dt2)

  d1_safe <- data.table::copy(dt1)
  d2_safe <- data.table::copy(dt2)

  # 2. Namespace Orthogonalization
  # Only suffix if specified and columns exist
  covs_1 <- setdiff(names(d1_safe), topo_keys)
  if (length(covs_1) > 0 && !is.null(suffix1)) {
    # Check if they are already suffixed to avoid d1_d1_d1
    to_suffix <- covs_1[!grepl(paste0("_", suffix1, "$"), covs_1)]
    if (length(to_suffix) > 0) data.table::setnames(d1_safe, to_suffix, paste0(to_suffix, "_", suffix1))
  }

  covs_2 <- setdiff(names(d2_safe), topo_keys)
  if (length(covs_2) > 0 && !is.null(suffix2)) {
    to_suffix <- covs_2[!grepl(paste0("_", suffix2, "$"), covs_2)]
    if (length(to_suffix) > 0) data.table::setnames(d2_safe, to_suffix, paste0(to_suffix, "_", suffix2))
  }

  # 3. Extract purely topological bounds
  b1 <- d1_safe[, .(pair_id = as.integer(pair_id), start = time, end = time_new)]
  b2 <- d2_safe[, .(pair_id = as.integer(pair_id), start = time, end = time_new)]
  stacked <- data.table::rbindlist(list(b1, b2))

  # 4. Execute C++ Splintering
  grid <- get_union_bounds(
    pair_id    = stacked$pair_id,
    time_start = as.numeric(stacked$start),
    time_end   = as.numeric(stacked$end)
  )
  data.table::setDT(grid)
  data.table::setnames(grid, c("pair_id", "grid_start", "grid_end"))
  grid[, pair_id := as.integer(pair_id)]

  if (nrow(grid) == 0) {
    return(data.table::data.table())
  }

  # 5. Parallel Projections
  # We join each dataset independently to the splintered grid.
  # This avoids complex scope shadowing issues in sequential non-equi joins.
  res1 <- d1_safe[grid, on = .(pair_id, time <= grid_start, time_new >= grid_end)]
  res2 <- d2_safe[grid, on = .(pair_id, time <= grid_start, time_new >= grid_end)]

  # 6. Coalesce & Merge
  # Coalesce topological columns (preference given to Dataset 1,
  # but if one is NA, we take the other)
  topo_to_coalesce <- intersect(topo_keys, names(res1))
  topo_to_coalesce <- intersect(topo_to_coalesce, names(res2))

  for (col in topo_to_coalesce) {
    # Non-equi joins set time/time_new to grid_start/grid_end values
    # We use fcoalesce to fill gaps in topology
    res1[, (col) := data.table::fcoalesce(get(col), res2[[col]])]
  }

  # Add Dataset 2 covariates to Dataset 1 result
  covs2_to_add <- setdiff(names(res2), names(res1))
  if (length(covs2_to_add) > 0) {
    res1[, (covs2_to_add) := res2[, covs2_to_add, with = FALSE]]
  }

  final_dt <- res1
  data.table::setorder(final_dt, pair_id, time)

  existing_topo <- intersect(topo_keys, names(final_dt))
  other_cols <- setdiff(names(final_dt), existing_topo)
  setcolorder(final_dt, c(existing_topo, other_cols))

  # Conditional NA handling
  if (impute_zero) {
    cov_cols <- setdiff(names(final_dt), topo_keys)
    for (col in cov_cols) set(final_dt, which(is.na(final_dt[[col]])), col, 0)
    if ("event" %in% names(final_dt)) set(final_dt, which(is.na(final_dt$event)), "event", 0)
  }

  if (omit_na) final_dt <- stats::na.omit(final_dt)
  return(final_dt)
}


union_dyadic_datasets_list <- function(dt_list,
                                       topo_keys = c("pair_id", "from", "to", "time", "time_new", "event", "status", "from_avail", "to_avail"),
                                       impute_zero = FALSE, omit_na = FALSE) {
  if (length(dt_list) == 0) {
    return(data.table())
  }
  if (length(dt_list) == 1) {
    return(dt_list[[1]])
  }

  # Use Reduce to apply union_dyadic_datasets in a nested fashion
  # We use the names of the list as suffixes if available
  nm <- names(dt_list)
  if (is.null(nm)) nm <- paste0("s", seq_along(dt_list))

  res <- dt_list[[1]]
  # Suffix the first one if not already suffixed
  covs <- setdiff(names(res), topo_keys)
  if (length(covs) > 0) {
    # Check if they are already suffixed to avoid d1_d1
    to_suffix <- covs[!grepl(paste0("_", nm[1], "$"), covs)]
    if (length(to_suffix) > 0) data.table::setnames(res, to_suffix, paste0(to_suffix, "_", nm[1]))
  }

  for (i in 2:length(dt_list)) {
    res <- union_dyadic_datasets(res, dt_list[[i]],
      topo_keys = topo_keys, suffix1 = NULL, suffix2 = nm[i],
      impute_zero = impute_zero, omit_na = omit_na
    )
  }

  return(res)
}

#' Preprocess Model Terms across Multiple Event Streams
#'
#' @param preprocessed Standard output from `formula_preprocess`.
#' @param n_nodes Number of nodes.
#' @param verbose Logical; if TRUE, print progress.
#' @param directed Logical; if TRUE, the model is directed.
#' @param simultaneous_interactions Logical; if TRUE, multiple interactions are allowed.
#' @param build_time Numeric; time at which to start building the dataset.
#' @param max_time Numeric; if positive, events after this time are excluded.
#'   Defaults to \code{-1.0} (no upper limit).
#' @param model_type Either "dem" or "rem".
#' @param impute_zero Logical; if TRUE, replace NAs in covariates with 0.
#' @param omit_na Logical; if TRUE, call na.omit() on the final table.
#'
#' @return A data.table containing the unified preprocessed data.
#' @keywords internal
preprocess_multi_stream <- function(preprocessed, n_nodes, verbose, directed, simultaneous_interactions,
                                    build_time = NULL, max_time = -1.0, model_type = "dem",
                                    impute_zero = TRUE, omit_na = TRUE) {
  build_time <- if (is.null(build_time)) 0 else build_time
  max_time <- if (is.null(max_time)) 0 else max_time

  stream_list <- preprocessed$stream_list
  # If everything is on the same stream (main events), we can use the fast path
  all_null <- all(sapply(stream_list, is.null))
  if (all_null) {
    if (model_type == "dem") {
      res <- preprocess(
        edgelist = as.matrix(preprocessed$events),
        terms = as.character(unlist(preprocessed$term_names)),
        data_list = preprocessed$data_list,
        transformations = as.character(unlist(preprocessed$transformation_list)),
        n_nodes = n_nodes,
        verbose = verbose,
        directed = directed,
        simultaneous_interactions = simultaneous_interactions,
        window_map = preprocessed$window_map,
        build_time = build_time,
        max_time = max_time
      )
    } else {
      res <- preprocess_rem(
        edgelist = as.matrix(preprocessed$events),
        terms = as.character(unlist(preprocessed$term_names)),
        data_list = preprocessed$data_list,
        transformations = as.character(unlist(preprocessed$transformation_list)),
        n_nodes = n_nodes,
        verbose = verbose,
        directed = directed,
        window_map = preprocessed$window_map,
        build_time = build_time,
        max_time = max_time
      )
    }
    res <- data.table::as.data.table(res)
    topo_names_s <- c("time_new", "time", "pair_id", "status", "event", "from", "to", "from_avail", "to_avail")
    cov_names_s <- preprocessed$coef_names
    data.table::setnames(res, c(topo_names_s, cov_names_s))

    # Ensure build_time filtering in fast path
    if (!is.null(build_time) && build_time > 0) {
      res <- res[res$time >= build_time]
    }

    return(res)
  }

  # Multi-stream path
  # 1. Identify unique streams
  unique_streams <- list(preprocessed$events)
  stream_assignments <- rep(1, length(stream_list))

  for (i in seq_along(stream_list)) {
    if (!is.null(stream_list[[i]])) {
      found <- FALSE
      if (length(unique_streams) >= 2) {
        for (j in 2:length(unique_streams)) {
          if (identical(stream_list[[i]], unique_streams[[j]])) {
            stream_assignments[i] <- j
            found <- TRUE
            break
          }
        }
      }
      if (!found) {
        unique_streams[[length(unique_streams) + 1]] <- stream_list[[i]]
        stream_assignments[i] <- length(unique_streams)
      }
    }
  }

  # 1.5 Security Checks and Data Validation
  for (s in seq_along(unique_streams)) {
    curr_stream <- as.matrix(unique_streams[[s]])
    s_label <- if (s == 1) "main event stream" else paste0("covariate stream ", s - 1)

    # Check dimensions
    if (ncol(curr_stream) < 4) {
      stop(sprintf("Multi-stream error: %s must have at least 4 columns (time, from, to, type).", s_label))
    }

    # Check node indices
    if (nrow(curr_stream) > 0) {
      nodes_in_stream <- unique(c(curr_stream[, 2], curr_stream[, 3]))
      if (any(nodes_in_stream > n_nodes) || any(nodes_in_stream < 1)) {
        stop(sprintf("Multi-stream error: %s contains node indices outside the range [1, %d].", s_label, n_nodes))
      }

      # For undirected events, we warn if the covariate stream seems to be directed but the model is undirected.
      # While preprocess() handles this by sorting, it's good practice to inform the user.
      if (!directed && any(curr_stream[, 2] > curr_stream[, 3])) {
        if (s > 1 && verbose) {
          # Only warn for covariate streams; the main stream is handled by dem()/rem()
          message(sprintf("Note: %s contains events (i, j) where i > j, but the model is undirected. These will be treated as undirected.", s_label))
        }
      }
    }
  }

  # 2. Preprocess each stream separately
  datasets <- list()
  for (s in seq_along(unique_streams)) {
    indices <- which(stream_assignments == s)
    # If s=1 (main events) and no terms assigned, we still need it for the topology/intervals
    # but we will have no covariate columns. C++ preprocess handles empty terms.

    terms_s <- preprocessed$term_names[indices]
    data_list_s <- preprocessed$data_list[indices]
    trans_s <- preprocessed$transformation_list[indices]

    # Inject a dummy type-3 event at build_time to ensure the state after the last
    # interaction is captured and extended appropriately for multi-stream splintering.
    edgelist_s <- as.matrix(unique_streams[[s]])
    if (!is.null(build_time) && nrow(edgelist_s) > 0) {
      max_t <- max(edgelist_s[, 1])
      if (max_t < build_time) {
        # Add type 3 event (exogenous change) at build_time.
        # from=1, to=2 are dummy nodes; type 3 triggers stat recording for all pairs in preprocess.
        dummy_event <- matrix(c(build_time, 1, 2, 3), nrow = 1)
        edgelist_s <- rbind(edgelist_s, dummy_event)
        edgelist_s <- edgelist_s[order(edgelist_s[, 1]), , drop = FALSE]
      }
    }


    if (model_type == "dem") {
      dt_s <- preprocess(
        edgelist = edgelist_s,
        terms = as.character(unlist(terms_s)),
        data_list = data_list_s,
        transformations = as.character(unlist(trans_s)),
        n_nodes = n_nodes,
        verbose = verbose,
        directed = directed,
        simultaneous_interactions = if (s == 1) simultaneous_interactions else TRUE,
        window_map = if (length(preprocessed$window_map) > 0) preprocessed$window_map else NULL,
        build_time = build_time,
        max_time = max_time
      )
    } else {
      dt_s <- preprocess_rem(
        edgelist = edgelist_s,
        terms = as.character(unlist(terms_s)),
        data_list = data_list_s,
        transformations = as.character(unlist(trans_s)),
        n_nodes = n_nodes,
        verbose = verbose,
        directed = directed,
        window_map = if (length(preprocessed$window_map) > 0) preprocessed$window_map else NULL,
        build_time = build_time,
        max_time = max_time
      )
    }

    # Handling empty returns (e.g., if a stream has no associated terms)
    if (length(dt_s) == 0 || nrow(dt_s) == 0) {
      if (length(indices) == 0) {
        # No terms for this stream - return an empty container that union_dyadic_datasets can ignore
        next
      }
      stop(sprintf("Stream %d returned no data but has %d terms. Interaction sequence might be invalid.", s, length(indices)))
    }

    dt_s <- data.table::as.data.table(dt_s)
    topo_names_s <- c("time_new", "time", "pair_id", "status", "event", "from", "to", "from_avail", "to_avail")
    cov_names_s <- preprocessed$coef_names[indices]
    data.table::setnames(dt_s, c(topo_names_s, cov_names_s))

    datasets[[s]] <- dt_s
  }

  # 3. Union all datasets
  names(datasets) <- paste0("s", seq_along(datasets))
  final_dt <- union_dyadic_datasets_list(datasets, impute_zero = impute_zero, omit_na = omit_na)

  # Explicitly filter by build_time for the multi-stream path
  if (!is.null(build_time) && build_time > 0) {
    if (nrow(final_dt) > 0) {
      final_dt <- final_dt[final_dt$time >= build_time]
    }
  }

  # 4. Filter and reorder columns
  all_coef_names <- names(preprocessed$coef_names)
  suffixed_names <- character(length(all_coef_names))
  for (i in seq_along(all_coef_names)) {
    s_id <- stream_assignments[i]
    suffixed_names[i] <- paste0(all_coef_names[i], "_s", s_id)
  }

  topo_keys_final <- c("time_new", "time", "pair_id", "status", "event", "from", "to", "from_avail", "to_avail")

  # Rename columns back to their original (unsuffixed) names for compatibility
  # but only for those that exist in final_dt
  existing_indices <- which(suffixed_names %in% names(final_dt))
  if (length(existing_indices) > 0) {
    data.table::setnames(final_dt, suffixed_names[existing_indices], all_coef_names[existing_indices])
  }

  final_cols <- c(topo_keys_final, all_coef_names[existing_indices])
  res <- final_dt[, final_cols, with = FALSE]
  return(res)
}


#' Internal helper to calculate predictions for gof
#'
#' @param model A redeem_result object.
#' @param data A data.table containing the covariates.
#' @keywords internal
calculate_predictions_helper <- function(model, data) {
  if (is.null(data) || nrow(data) == 0) {
    return(numeric(0))
  }

  # coefficients names from the model
  nm <- names(model$coefficients)

  # Identify columns
  is_degree <- grep("^effect_|^sender_|^receiver_", nm)
  is_time <- grep("^time_cat", nm)
  is_core <- setdiff(seq_along(nm), c(is_degree, is_time))

  # Reconstruct linear predictor
  # 1. Core (relational) effects
  pred <- rep(0, nrow(data))
  if (length(is_core) > 0) {
    # Check which columns exist in data
    existing_cols <- intersect(nm[is_core], names(data))
    if (length(existing_cols) > 0) {
      # Use matrix multiplication for speed
      X <- as.matrix(data[, existing_cols, with = FALSE])
      beta <- model$coefficients[existing_cols]
      pred <- as.vector(X %*% beta)
    }
  }

  # 2. Add degree effects
  # 2. Add degree effects
  # Use est_degree directly rather than coefficient-name matching, because
  # standardized redeem results may strip degree terms from coefficients.
  deg_vals <- model$est_degree
  if (!is.null(deg_vals) && length(deg_vals) > 0) {
    # In data, we have 'from' and 'to' columns.
    observed_n_nodes <- max(c(data$from, data$to), na.rm = TRUE)

    # Prefer stored metadata when available, otherwise infer it.
    directed <- model$directed
    if (is.null(directed)) {
      directed <- length(deg_vals) >= 2 * observed_n_nodes
    }

    n_nodes <- model$n_nodes
    if (is.null(n_nodes)) {
      if (directed) {
        n_nodes <- min(observed_n_nodes, floor(length(deg_vals) / 2))
      } else {
        n_nodes <- min(observed_n_nodes, length(deg_vals))
      }
    }

    if (directed) {
      # 1:n_nodes is mu (sender), n_nodes+1:2n_nodes is nu (receiver)
      if (length(deg_vals) >= 2 * n_nodes) {
        mu <- deg_vals[1:n_nodes]
        nu <- deg_vals[(n_nodes + 1):(2 * n_nodes)]
        pred <- pred + mu[data$from] + nu[data$to]
      }
    } else {
      if (length(deg_vals) >= n_nodes) {
        pred <- pred + deg_vals[data$from] + deg_vals[data$to]
      }
    }
  }

  # 3. Add temporal effects
  est_time <- model$est_time
  time_changepoints <- model$time_changepoints

  if (!is.null(est_time) && length(est_time) > 0 && !is.null(time_changepoints)) {
    # Blockwise/MM path: Use integrated intensity via get_time_offset
    # Reconstruct topological slice indices for the integration helper
    time_slices <- findInterval(data$time_new, c(-Inf, time_changepoints, Inf))
    time_slices_from <- findInterval(data$time, c(-Inf, time_changepoints, Inf))

    # Temporal multiplier per row (already accounted for baseline bin = 0)
    baseline_mult <- as.vector(get_time_offset(
      from_slice_r = time_slices_from,
      to_slice_r = time_slices,
      from_time_r = data$time,
      to_time_r = data$time_new,
      est_time_r = if (isTRUE(model$full_baseline)) exp(as.vector(est_time)) else exp(c(0, as.vector(est_time))),
      changepoints_r = c(time_changepoints, max(data$time_new, na.rm = TRUE) + 1e-6)
    ))

    # Intensity: exp(pred) * baseline_mult
    res <- exp(pred) * baseline_mult
  } else if (length(is_time) > 0) {
    # Parametric/NR path: time_cat columns are indicators (1/0)
    existing_time <- intersect(nm[is_time], names(data))
    if (length(existing_time) > 0) {
      X_time <- as.matrix(data[, existing_time, with = FALSE])
      beta_time <- model$coefficients[existing_time]
      pred <- pred + as.vector(X_time %*% beta_time)
    }
    # Final intensity: exp(pred + offset)
    res <- exp(pred + data$offset)
  } else {
    # No temporal effects beyond baseline offset
    res <- exp(pred + data$offset)
  }

  return(res)
}

#' Reconstruct estimation data from a model object
#'
#' @param models A redeem model object (rem or dem).
#' @param verbose Logical.
#' @keywords internal
reproduce_model_data <- function(models, verbose = FALSE) {
  # 1. Collect metadata
  model_0_1 <- if (inherits(models, "rem")) models$model else models$model_0_1
  model_1_0 <- if (inherits(models, "rem")) NULL else models$model_1_0

  formula_0_1 <- if (inherits(models, "rem")) models$formula else models$formula_0_1
  formula_1_0 <- if (inherits(models, "rem")) NULL else models$formula_1_0

  # n_nodes and other globals
  n_nodes <- models$n_nodes
  directed <- models$directed
  build_time <- if (is.null(models$build_time)) 0 else models$build_time
  max_time <- if (is.null(models$max_time)) 0 else models$max_time

  # 2. Re-preprocess
  preprocessed <- if (!is.null(models$preprocessed)) {
    models$preprocessed
  } else {
    formula_preprocess(
      formula_1_0 = if (inherits(models, "rem")) models$formula else formula_1_0,
      formula_0_1 = if (inherits(models, "rem")) models$formula else formula_0_1,
      events = models$events,
      n_nodes = n_nodes,
      model_type = if (inherits(models, "rem")) "rem" else "dem",
      directed = directed
    )
  }

  # 3. Reconstruct full splintered timeline
  data_reproduced <- preprocess_multi_stream(
    preprocessed = preprocessed,
    model_type = if (inherits(models, "rem")) "rem" else "dem",
    n_nodes = n_nodes,
    directed = directed,
    simultaneous_interactions = if (is.null(models$simultaneous_interactions)) FALSE else models$simultaneous_interactions,
    build_time = build_time,
    max_time = max_time,
    verbose = verbose
  )

  # 4. Cleanup and offset (Mirroring rem.R/dem.R)
  data.table::setDT(data_reproduced)
  data_reproduced[, diff := time_new - time]
  data_reproduced[diff <= 0, diff := NA]
  data_reproduced[, offset := log(diff)]
  data_reproduced[is.infinite(offset), event := NA]
  data_reproduced <- data_reproduced[!is.na(event)]

  # 5. Temporal binning
  # Use changepoints from the newly generated preprocessed object.
  # For DEM objects, formation and dissolution may have different baseline
  # changepoints/labels, so construct time bins separately by process.
  add_time_binning <- function(dt, time_changepoints, labels_changepoints) {
    if (is.null(time_changepoints) || length(time_changepoints) == 0) {
      return(dt)
    }
    if (is.null(labels_changepoints)) {
      labels_changepoints <- as.character(time_changepoints)
    }
    full_labels <- c("Beg", labels_changepoints)
    dt[, time_cat := cut(time_new,
      breaks = c(-1, time_changepoints, Inf),
      labels = full_labels
    )]
    time_mm <- model.matrix(~time_cat, dt)
    cbind(dt, time_mm)
  }

  if (inherits(models, "rem")) {
    data_reproduced <- add_time_binning(
      data_reproduced,
      preprocessed$baseline_changepoints_0_1,
      preprocessed$baseline_labels_0_1
    )
  } else {
    data_0_1_tmp <- data_reproduced[status == 0]
    data_1_0_tmp <- data_reproduced[status == 1]

    data_0_1_tmp <- add_time_binning(
      data_0_1_tmp,
      preprocessed$baseline_changepoints_0_1,
      preprocessed$baseline_labels_0_1
    )
    data_1_0_tmp <- add_time_binning(
      data_1_0_tmp,
      preprocessed$baseline_changepoints_1_0,
      preprocessed$baseline_labels_1_0
    )

    data_reproduced <- data.table::rbindlist(
      list(data_0_1_tmp, data_1_0_tmp),
      use.names = TRUE,
      fill = TRUE
    )
  }

  # 6. Add predictions back
  if (inherits(models, "rem")) {
    data_reproduced$prediction <- calculate_predictions_helper(models$model, data_reproduced)
    return(data_reproduced)
  } else {
    # For DEM, split by process
    data_0_1 <- data_reproduced[status == 0]
    data_1_0 <- data_reproduced[status == 1]

    if (!is.null(models$model_0_1)) {
      data_0_1$prediction <- calculate_predictions_helper(models$model_0_1, data_0_1)
    }
    if (!is.null(models$model_1_0)) {
      data_1_0$prediction <- calculate_predictions_helper(models$model_1_0, data_1_0)
    }

    return(list(data_0_1 = data_0_1, data_1_0 = data_1_0))
  }
}

#' Match user-provided coefficients to internal model names
#'
#' @param user_coefs Named or unnamed vector of coefficients provided by the user.
#' @param internal_names Vector of internal coefficient names (labels).
#' @param internal_keys Vector of internal coefficient keys (e.g., intercept, dyadic_cov_identity).
#'
#' @return A numeric vector of the same length as internal_names.
#' @keywords internal
match_coefficients <- function(user_coefs, internal_names, internal_keys = NULL) {
  n_internal <- length(internal_names)
  n_user <- length(user_coefs)
  res <- numeric(n_internal)

  if (n_user == 0) {
    return(res)
  }

  user_names <- names(user_coefs)

  # If no names provided at all, use positional matching
  if (is.null(user_names)) {
    if (n_user == n_internal) {
      return(as.numeric(user_coefs))
    }
    # If user provides k-1 and internal has Intercept at start, offset by 1
    if (n_user == (n_internal - 1) && !is.null(internal_names) && tolower(internal_names[1]) == "intercept") {
      res[2:n_internal] <- as.numeric(user_coefs)
      return(res)
    }
    # Otherwise, fill from start
    n_fill <- min(n_user, n_internal)
    res[1:n_fill] <- as.numeric(user_coefs[1:n_fill])
    return(res)
  }

  # Named matching
  # 1. Try matching against internal keys first (e.g. inertia_log)
  idx <- rep(NA, n_user)
  if (!is.null(internal_keys)) {
    idx <- match(user_names, internal_keys)
  }

  # 2. Try matching against internal names (labels like "Intercept")
  unmatched <- which(is.na(idx) & user_names != "")
  if (length(unmatched) > 0) {
    # Ensure internal_names is character vector before match
    internal_names_char <- as.character(internal_names)
    m <- match(user_names[unmatched], internal_names_char)
    idx[unmatched] <- m
  }

  # 3. Fuzzy match for still unmatched named elements (ignoring suffixes)
  unmatched <- which(is.na(idx) & user_names != "")
  if (length(unmatched) > 0) {
    clean_internal <- gsub("_log$|_identity$|_bin$|_sig$|_recip$", "", internal_keys)
    if (is.null(clean_internal) || length(clean_internal) == 0) {
      clean_internal <- gsub("\\s*\\(transformation = .*\\)$", "", internal_names)
      clean_internal <- gsub("\\s*\\(.*\\)$", "", clean_internal)
    }

    # Handle plural/singular differences (e.g. current_common_partners vs current_common_partner)
    clean_internal_plural <- gsub("s$", "", clean_internal)
    user_names_plural <- gsub("s$", "", user_names)

    for (i in unmatched) {
      m <- match(user_names[i], clean_internal)
      if (is.na(m)) m <- match(user_names_plural[i], clean_internal)
      if (is.na(m)) m <- match(user_names[i], clean_internal_plural)
      if (is.na(m)) m <- match(user_names_plural[i], clean_internal_plural)

      if (is.na(m)) {
        # Check if user name is contained in any internal name (case insensitive)
        contains <- grepl(user_names[i], internal_names, ignore.case = TRUE)
        if (any(contains)) m <- which(contains)[1]
      }

      if (!is.na(m)) idx[i] <- m
    }
  }

  # 4. Assign matched elements
  matched <- !is.na(idx)
  assigned_internal_slots <- numeric(0)
  if (any(matched)) {
    valid_slots <- idx[matched]
    res[valid_slots] <- as.numeric(user_coefs[matched])
    assigned_internal_slots <- valid_slots
  }

  # 5. Handle unnamed or non-matched elements positionally ONLY IF no names were provided at all
  # OR if we explicitly want to allow it (not recommended for named input)
  remaining_user_idx <- which(!matched)
  if (length(remaining_user_idx) > 0 && is.null(user_names)) {
    remaining_internal_idx <- setdiff(1:n_internal, assigned_internal_slots)
    n_fill <- min(length(remaining_user_idx), length(remaining_internal_idx))
    if (n_fill > 0) {
      res[remaining_internal_idx[1:n_fill]] <- as.numeric(user_coefs[remaining_user_idx[1:n_fill]])
    }
  }

  return(res)
}

#' Control Parameters for REDEEM Models
#'
#' Unified control object to manage estimation parameters for \code{\link{rem}}
#' and \code{\link{dem}} functions.
#'
#' @param it_max Integer; maximum number of iterations for the algorithm.
#'   Defaults to 100.
#' @param tol Numeric; convergence tolerance. Defaults to 1e-10.
#' @param accelerated Logical; if \code{TRUE}, uses SQUAREM acceleration for
#'   MM updates. Defaults to FALSE.
#' @param verbose Logical; if \code{TRUE}, prints progress information. Defaults to FALSE.
#' @param weighting Logical; whether to use weighting to group identical
#'   observations. Defaults to TRUE.
#' @param subsample Numeric; proportion of data to subsample for internal
#'   GLM checks. Defaults to 1.
#' @param use_glm Logical; if \code{TRUE}, uses standard GLM for updating
#'   core coefficients. This is often slower but can yield more robust updates.
#'   Defaults to FALSE.
#' @param legacy Logical; if \code{TRUE}, uses
#'   a single \code{glm.fit} call instead of the iterative loop.
#'   Defaults to FALSE.
#' @param build_time Numeric; time at which to start building the estimation
#'   dataset. Events before this time are used to compute statistics but not
#'   included as observations. Defaults to NULL, in which case all events are included.
#' @param return_data Logical; whether to return preprocessed data frames
#'   in the result. Defaults to FALSE.
#' @param save_hist Logical; whether to save the iteration history of
#'   coefficients. Defaults to TRUE.
#' @param estimate Character; estimation method for \code{\link{dem}} and
#'   \code{\link{rem}} ("Blockwise", "NR", or "GD"). Defaults to "Blockwise".
#' @param check_matrix Logical; whether to apply \code{\link{check_matrix}}
#'   to the event data before estimation. If \code{TRUE}, repairs missing
#'   events (e.g., adding start events for interactions that only have end
#'   events). Defaults to FALSE.
#' @param inf_unidentifiable Logical; whether to set unidentifiable
#'   coefficients (e.g., actors with 0 event counts, globally
#'   invariant/collinear covariates) to \code{-Inf}. Defaults to TRUE.
#'
#' @return A list of class \code{"redeem_control"} containing the specified
#'   parameters.
#'
#' @export
control.redeem <- function(it_max = 100,
                           tol = 1e-10,
                           accelerated = FALSE,
                           verbose = FALSE,
                           weighting = TRUE,
                           subsample = 1,
                           build_time = NULL,
                           # use_cpp = TRUE,
                           use_glm = FALSE,
                           return_data = FALSE,
                           save_hist = TRUE,
                           estimate = "Blockwise",
                           legacy = FALSE,
                           check_matrix = FALSE,
                           inf_unidentifiable = TRUE) {

  # Parameter validation
  if (!is.numeric(it_max) || any(it_max <= 0) || length(it_max) > 2) {
    stop("it_max must be a positive integer (length 1 or 2).")
  }
  if (!is.numeric(tol) || any(tol <= 0) || length(tol) > 2) {
    stop("tol must be a positive numeric value (length 1 or 2).")
  }
  if (!is.numeric(subsample) || subsample <= 0 || subsample > 1) {
    stop("subsample must be between 0 and 1.")
  }
  if (!is.logical(accelerated) || length(accelerated) > 2) {
    stop("accelerated must be logical (length 1 or 2).")
  }
  if (!is.logical(verbose)) stop("verbose must be logical.")
  if (!is.logical(weighting)) stop("weighting must be logical.")
  # if (!is.logical(use_cpp)) stop("use_cpp must be logical.")
  if (!is.logical(use_glm)) stop("use_glm must be logical.")
  if (!is.logical(return_data)) stop("return_data must be logical.")
  if (!is.logical(save_hist)) stop("save_hist must be logical.")
  if (!is.logical(check_matrix)) stop("check_matrix must be logical.")
  if (!is.logical(legacy)) stop("legacy must be logical.")
  if (!is.logical(inf_unidentifiable)) stop("inf_unidentifiable must be logical.")


  if (!estimate %in% c("Blockwise", "NR", "GD")) {
    stop("estimate must be either 'Blockwise', 'NR', or 'GD'.")
  }

  res <- list(
    it_max = as.integer(it_max),
    tol = as.numeric(tol),
    accelerated = accelerated,
    verbose = verbose,
    weighting = weighting,
    subsample = subsample,
    use_cpp = TRUE,
    use_glm = use_glm,
    return_data = return_data,
    save_hist = save_hist,
    estimate = estimate,
    build_time = build_time,
    legacy = legacy,
    check_matrix = check_matrix,
    inf_unidentifiable = inf_unidentifiable
  )

  class(res) <- "redeem_control"
  return(res)
}

#' Core Estimation Logic for REM and DEM Transitions
#'
#' This internal helper function encapsulates the estimation routines for a single
#' relational or durational event transition. It is used by both `rem()` and `dem()`.
#'
#' @import Rcpp
#' @keywords internal
estimate_transition <- function(data,
                                formula_original,
                                formula_new,
                                indicators,
                                n_nodes,
                                estimate_method,
                                it_max,
                                tol,
                                accelerated,
                                subsample,
                                verbose,
                                estimate_degree,
                                directed,
                                semiparametric = FALSE,
                                labels_changepoints = NULL,
                                time_changepoints = NULL,
                                coef_init = NULL,
                                model_type = "dem",
                                process = "0-1",
                                return_data = TRUE, save_hist = TRUE,
                                use_glm = FALSE,
                                legacy = FALSE,
                                inf_unidentifiable = TRUE,
                                events = NULL) {
  if (semiparametric) {
    # Proportional intensity estimation via Cox model
    model <- survival::coxph(
      update(formula_new, new = "survival::Surv(time = time, time2 = time_new, event = event) ~ . - Intercept"),
      data = data
    )
    return(model)
  }

  # Initial coefficients
  if (is.null(coef_init)) {
    start_coefs <- get_start_coefs(
      formula = formula_original,
      n_nodes = n_nodes,
      fixed_effects = estimate_degree,
      time_changepoints = time_changepoints,
      directed = directed,
      model_type = model_type,
      process = process
    )
    coefs_core <- start_coefs$coef_core
    coefs_degree <- start_coefs$coef_degree
    coefs_time <- start_coefs$coef_time
  } else {
    # Parse provided coefficients robustly
    start_info <- get_start_coefs(
      formula = formula_original,
      n_nodes = n_nodes,
      fixed_effects = estimate_degree,
      time_changepoints = time_changepoints,
      directed = directed,
      model_type = model_type,
      process = process
    )
    number_suff <- length(start_info$coef_core)
    number_time <- length(start_info$coef_time)
    expected_degree <- if (directed) 2 * n_nodes else n_nodes
    if (length(coef_init) != (number_suff + expected_degree + number_time)) {
      # Fallback to get_start_coefs if size mismatch, but usually we stop
      stop(paste0(
        "Wrong number of coefficients provided for initialization. Expected ",
        (number_suff + expected_degree + number_time), " but got ", length(coef_init)
      ))
    }
    coefs_core <- coef_init[seq_len(number_suff)]
    coefs_degree <- coef_init[(number_suff + 1):(number_suff + expected_degree)]
    coefs_time <- if (number_time > 0) {
      coef_init[(number_suff + expected_degree + 1):(number_suff + expected_degree + number_time)]
    } else {
      numeric(0)
    }
  }

  if (estimate_method == "NR") {
    # Full Newton-Raphson using glm.fit
    formula_tmp <- update(formula_new, new = "event ~ . -1")
    environment(formula_tmp) <- environment()
    mm <- tryCatch(
      {
        model.matrix(formula_tmp, data)
      },
      error = function(e) {
        stop(paste("Error creating model matrix:", e$message))
      }
    )

    # Newton-Raphson assembly: combine all relevant coefficients
    coef_current <- c(
      if (!is.null(coefs_core)) coefs_core else numeric(0),
      if (estimate_degree && !is.null(coefs_degree)) coefs_degree else numeric(0),
      if (!is.null(coefs_time)) coefs_time else numeric(0)
    )
    if (length(coef_current) == 0) coef_current <- NULL

    if (legacy) {
      cat("Legacy")
      # Perfectly match DEM package: single glm call
      model <- tryCatch(
        {
          suppressWarnings(stats::glm(
            formula = formula_tmp,
            data = data,
            family = stats::poisson(),
            offset = data$offset,
            weights = data$weight,
            start = coef_current,
            control = stats::glm.control(maxit = it_max, epsilon = tol)
          ))
        },
        error = function(e) {
          suppressWarnings(stats::glm(
            formula = formula_tmp,
            data = data,
            family = stats::poisson(),
            offset = data$offset,
            weights = data$weight,
            control = stats::glm.control(maxit = it_max, epsilon = tol)
          ))
        }
      )
      it <- if (!is.null(model$iter)) model$iter else 1
      coef_hist <- matrix(model$coefficients, nrow = 1)
      colnames(coef_hist) <- names(model$coefficients)
      intensity_nr <- model$fitted.values / exp(if (!is.null(data$offset)) data$offset else 0)
      llh_hist <- calc_llh_scaled(
        pred = model$fitted.values,
        intensity = intensity_nr,
        delta = data$event,
        pair_id = if (!is.null(data$pair_id)) data$pair_id else seq_len(nrow(data))
      )
    } else {
      # Full Newton-Raphson using manual loop for better output control
      dev_old <- Inf
      llh_hist <- numeric(it_max)
      coef_hist <- matrix(NA, nrow = it_max, ncol = if (length(coef_current) == 0) ncol(mm) else length(coef_current))

      for (it in 1:it_max) {
        model_res <- tryCatch(
          {
            suppressWarnings(stats::glm.fit(
              x = mm,
              y = data$event,
              family = stats::poisson(),
              offset = data$offset,
              weights = data$weight,
              start = coef_current,
              control = stats::glm.control(maxit = 1)
            ))
          },
          error = function(e) {
            suppressWarnings(stats::glm.fit(
              x = mm,
              y = data$event,
              family = stats::poisson(),
              offset = data$offset,
              weights = data$weight,
              control = stats::glm.control(maxit = 1)
            ))
          }
        )

        if (is.null(model_res) || !is.list(model_res)) {
          stop("Failed to estimate model using NR (glm.fit failed). This often happens with singular matrices or invalid data.")
        }
        model <- model_res

        coef_current <- model$coefficients
        coef_hist[it, ] <- coef_current
        intensity_nr <- model$fitted.values / exp(if (!is.null(data$offset)) data$offset else 0)
        llh_hist[it] <- calc_llh_scaled(
          pred = model$fitted.values,
          intensity = intensity_nr,
          delta = data$event,
          pair_id = if (!is.null(data$pair_id)) data$pair_id else seq_len(nrow(data))
        )


        if (verbose) {
          val_check <- if (it > 1) abs(model$deviance - dev_old) else Inf
          cat("\rIteration:", it, "/", it_max, "- LLH:", round(llh_hist[it], 4), "- Criterion:", round(val_check, 4))
          utils::flush.console()
        }
        if (it > 1 && abs(model$deviance - dev_old) < 1e-8) {
          if (verbose) cat("\nNR Converged after", it, "iterations\n")
          break
        }
        dev_old <- model$deviance
      }
      llh_hist <- llh_hist[1:it]
      coef_hist <- coef_hist[1:it, , drop = FALSE]
    }

    idx_pos <- data$event > 0
    is_separated <- colSums(mm[idx_pos, , drop = FALSE] != 0) == 0
    if (inf_unidentifiable) {
      model$coefficients[is_separated] <- -Inf
    }

    # Standardize the output using the common helper
    nm <- names(model$coefficients)
    is_degree <- grep("^effect_|^sender_|^receiver_", nm)
    is_time <- grep("^time_cat", nm)
    is_core <- setdiff(seq_along(nm), c(is_degree, is_time))

    # Coefficients for the main table (fixed effects + temporal bins)
    idx_main <- c(is_core, is_time)

    # Handle covariance for NR
    covariance_main <- if (!is.null(model$qr)) {
      p <- model$rank
      p1 <- 1L:p
      Qr <- model$qr
      tryCatch(
        {
          full_cov <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
          colnames(full_cov) <- rownames(full_cov) <- nm[p1]
          idx_cov <- intersect(nm[p1], nm[idx_main])
          full_cov[idx_cov, idx_cov, drop = FALSE]
        },
        error = function(e) matrix(NA, length(idx_main), length(idx_main))
      )
    } else {
      matrix(NA, length(idx_main), length(idx_main))
    }

    # Set covariance to NA for separated terms
    sep_names <- colnames(mm)[is_separated]
    if (length(sep_names) > 0) {
      idx_sep_cov <- intersect(colnames(covariance_main), sep_names)
      if (length(idx_sep_cov) > 0) {
        covariance_main[idx_sep_cov, ] <- NA
        covariance_main[, idx_sep_cov] <- NA
      }
    }

    # Assign column names to history matrix for plotting
    col_names_final <- colnames(mm)
    colnames(coef_hist) <- col_names_final

    has_intercept <- any(tolower(col_names_final) %in% c("intercept", "(intercept)"))
    full_baseline <- !(estimate_degree || has_intercept)

    res <- standardize_redeem_result(
      coefficients = model$coefficients,
      coef_hist = coef_hist,
      covariance = covariance_main,
      llh_hist = llh_hist,
      data = data,
      prediction = model$fitted.values,
      method = "dem.nr",
      n_nodes = n_nodes,
      directed = directed,
      time_changepoints = time_changepoints,
      labels_changepoints = labels_changepoints,
      full_baseline = full_baseline,
      return_data = return_data, save_hist = save_hist
    )
    return(res)
  } else if (estimate_method == "GD") {
    # Gradient Descent (not optimized currently)
    X <- model.matrix(update(formula_new, new = "event ~ . + offset(offset) - 1"), data = data)

    # Ensure coef matches X columns
    all_coefs <- c(coefs_core, coefs_degree, coefs_time)
    # If the lengths don't match, we need to subset or use start_coefs logic
    if (length(all_coefs) != ncol(X)) {
      # Check if we can match by name or just use ncol(X)
      coef_input <- rep(0, ncol(X))
      names(coef_input) <- colnames(X)
      # Try to fill in from all_coefs if named, else just keep zeros
      # For now, simplest is to just use a zero vector of correct length if mismatch
      # as GD is sensitive to starting points anyway
    } else {
      coef_input <- all_coefs
    }

    model <- gd_estimation(
      X = X,
      y = data$event,
      coef = coef_input,
      offset = data$offset,
      max_iter = it_max,
      tol = tol
    )

    # Identify separated coefficients (MLE -> -Inf)
    idx_pos <- data$event > 0
    is_separated <- colSums(X[idx_pos, , drop = FALSE] != 0) == 0
    if (inf_unidentifiable) {
      model$coef[is_separated] <- -Inf
    }

    # Standardize output using the common helper
    nm <- colnames(X)
    names(model$coef) <- nm

    # Calculate covariance for non-pop effects at the end of GD
    pred <- as.vector(exp(X %*% model$coef + data$offset))
    is_degree <- grep("^effect_|^sender_|^receiver_", nm)
    is_time <- grep("^time_cat", nm)
    is_core <- setdiff(seq_along(nm), c(is_degree, is_time))

    W <- pred
    H <- t(X) %*% (W * X)
    full_cov <- tryCatch(
      {
        solve(H)
      },
      error = function(e) MASS::ginv(H)
    )
    colnames(full_cov) <- rownames(full_cov) <- nm
    idx_main <- c(is_core, is_time)
    covariance_main <- full_cov[idx_main, idx_main, drop = FALSE]

    # Set covariance to NA for separated terms
    sep_names <- colnames(X)[is_separated]
    if (length(sep_names) > 0) {
      idx_sep_cov <- intersect(colnames(covariance_main), sep_names)
      if (length(idx_sep_cov) > 0) {
        covariance_main[idx_sep_cov, ] <- NA
        covariance_main[, idx_sep_cov] <- NA
      }
    }

    # Assign column names to history matrix for plotting
    colnames(model$coef_hist) <- nm

    has_intercept <- any(tolower(nm) == "intercept")
    full_baseline <- !(estimate_degree || has_intercept)

    res <- standardize_redeem_result(
      coefficients = model$coef,
      coef_hist = model$coef_hist,
      covariance = covariance_main,
      llh_hist = model$llh,
      data = data,
      prediction = pred,
      method = "dem.gd",
      n_nodes = n_nodes,
      directed = directed,
      time_changepoints = time_changepoints,
      labels_changepoints = labels_changepoints,
      full_baseline = full_baseline,
      return_data = return_data, save_hist = save_hist
    )
    return(res)
  } else {
    # Blockwise MM/NR estimation
    if (is.null(time_changepoints) && length(time_changepoints) == 0) {
      est_time = NULL
      time_changepoints = NULL
      labels_changepoints = NULL
    }
    model <- estimate_mmt(
      data = data,
      indicators = indicators,
      it_max = it_max,
      n_nodes = n_nodes,
      tol = tol,
      accelerated = accelerated,
      labels_changepoints = labels_changepoints,
      time_changepoints = time_changepoints,
      subsample = subsample,
      verbose = verbose,
      est_degree = coefs_degree,
      est_core = coefs_core,
      est_time = coefs_time,
      estimate_degree = estimate_degree,
      directed = directed,
      return_data = return_data, save_hist = save_hist,
      use_glm = use_glm,
      inf_unidentifiable = inf_unidentifiable
    )
    return(model)
  }
}

#' Standardize Estimation Output for Redeem Models
#'
#' @param coefficients Joint coefficient vector.
#' @param coef_hist Matrix of coefficient history.
#' @param covariance Covariance matrix for core effects.
#' @param llh_hist Vector of log-likelihood history.
#' @param data Input data.
#' @param prediction Predicted values.
#' @param method Estimation method name (for class).
#' @param n_nodes Number of nodes.
#' @param directed Logical; are the interaction events directed?
#' @param time_changepoints Numeric vector of time changepoints.
#' @param labels_changepoints Character vector of labels for the changepoints.
#' @param return_data Logical; should the estimation dataset be returned?
#' @param save_hist Logical; should the parameter history be returned?
#' @keywords internal
standardize_redeem_result <- function(coefficients,
                                      coef_hist,
                                      covariance,
                                      llh_hist,
                                      data,
                                      prediction,
                                      method = "dem.nr",
                                      n_nodes = NULL,
                                      directed = FALSE,
                                      time_changepoints = NULL,
                                      labels_changepoints = NULL,
                                      full_baseline = FALSE,
                                      return_data = TRUE, save_hist = TRUE) {
  nm <- names(coefficients)
  if (is.null(nm)) nm <- colnames(coef_hist)

  is_degree <- grep("^effect_[0-9]+$|^sender_[0-9]+$|^receiver_[0-9]+$", nm)
  is_time <- grep("^time_cat", nm)
  is_core <- setdiff(seq_along(nm), c(is_degree, is_time))

  est_core <- coefficients[is_core]
  est_degree <- coefficients[is_degree]
  est_time <- coefficients[is_time]
  if (length(est_time) > 0) {
    names(est_time) <- gsub("^time_cat", "time_", names(est_time))
    names(coefficients)[is_time] <- names(est_time)
  }

  if (!is.null(covariance)) {
    curr_cov_nm <- colnames(covariance)
    if (!is.null(curr_cov_nm)) {
      colnames(covariance) <- rownames(covariance) <- gsub("^time_cat", "time_", curr_cov_nm)
    }
  }

  # Ensure data has prediction
  data$prediction <- prediction

  res <- list(
    coefficients = coefficients,
    est_degree = est_degree,
    est_degree_hist = if (save_hist && length(is_degree) > 0 && !is.null(coef_hist)) coef_hist[, is_degree, drop = FALSE] else NULL,
    est_core = est_core,
    coefficients_core_hist = if (save_hist && !is.null(coef_hist)) coef_hist[, is_core, drop = FALSE] else NULL,
    est_time = est_time,
    est_time_hist = if (save_hist && length(is_time) > 0 && !is.null(coef_hist)) {
      tmp_hist <- coef_hist[, is_time, drop = FALSE]
      colnames(tmp_hist) <- names(est_time)
      tmp_hist
    } else NULL,
    covariance = covariance,
    llh = llh_hist[length(llh_hist)],
    llh_hist = llh_hist,
    data = if (return_data) data else NULL,
    n_nodes = n_nodes,
    directed = directed,
    time_changepoints = time_changepoints,
    labels_changepoints = labels_changepoints,
    full_baseline = full_baseline,
    n_obs = nrow(data),
    prediction = prediction
  )
  class(res) <- c(method, "redeem_result")
  return(res)
}

#' Process event actor columns and automatically identify or validate n_nodes
#'
#' @param events A matrix or data frame of events with columns \code{from} and \code{to} (or columns 2 and 3).
#' @param n_nodes Integer; the total number of actors in the network, or \code{NULL}.
#' @param directed Logical; whether the interaction events are directed. Defaults to TRUE.
#'
#' @return A list containing \code{events} (potentially modified) and \code{n_nodes}.
#' @keywords internal
process_event_actors <- function(events, n_nodes = NULL, directed = TRUE) {
  col_from <- if ("from" %in% colnames(events)) "from" else 2
  col_to <- if ("to" %in% colnames(events)) "to" else 3

  if (is.data.frame(events)) {
    from_vals <- events[[col_from]]
    to_vals <- events[[col_to]]
  } else {
    from_vals <- events[, col_from]
    to_vals <- events[, col_to]
  }

  is_string <- is.character(from_vals) || is.factor(from_vals) || 
               is.character(to_vals) || is.factor(to_vals)

  if (is_string) {
    from_vals <- as.character(from_vals)
    to_vals <- as.character(to_vals)
  }

  unique_actors <- unique(c(from_vals, to_vals))
  unique_actors <- unique_actors[!is.na(unique_actors) & unique_actors != "NA"]
  K <- length(unique_actors)

  if (!is.null(n_nodes)) {
    if (is_string) {
      if (K > n_nodes) {
        stop(sprintf("More unique actors (%d) identified in the events than provided in the n_nodes argument (%d).", K, n_nodes))
      }
    } else {
      # For numeric actors:
      max_actor_id <- as.integer(max(c(from_vals, to_vals), na.rm = TRUE))
      if (K > n_nodes) {
        stop(sprintf("More unique actors (%d) identified in the events than provided in the n_nodes argument (%d).", K, n_nodes))
      }
      if (max_actor_id > n_nodes) {
        stop(sprintf("Maximum actor ID (%d) is larger than the provided n_nodes argument (%d).", max_actor_id, n_nodes))
      }
    }
  }

  if (is_string) {
    unique_actors_sorted <- sort(unique_actors)
    from_vals <- match(from_vals, unique_actors_sorted)
    to_vals <- match(to_vals, unique_actors_sorted)
    if (is.null(n_nodes)) {
      n_nodes <- K
    }
  } else {
    if (is.null(n_nodes)) {
      n_nodes <- as.integer(max(c(from_vals, to_vals), na.rm = TRUE))
    }
  }

  if (!directed) {
    swap <- from_vals > to_vals
    if (any(swap, na.rm = TRUE)) {
      tmp <- from_vals[swap]
      from_vals[swap] <- to_vals[swap]
      to_vals[swap] <- tmp
    }
  }

  if (is.data.frame(events)) {
    events[[col_from]] <- from_vals
    events[[col_to]] <- to_vals
  } else {
    events[, col_from] <- from_vals
    events[, col_to] <- to_vals
    if (is.matrix(events) || is.array(events)) {
      storage.mode(events) <- "numeric"
    }
  }

  return(list(events = events, n_nodes = n_nodes))
}

