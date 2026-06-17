#' Initialize redeem Model Terms
#'
#' This is an internal dispatcher that calls the appropriate term
#' initialization function.
#'
#' @param term_name The name of the term.
#' @param arglist A list of arguments passed to the term in the formula.
#' @param model_type Either "dem" or "rem".
#' @param process Either "0-1" (incidence) or "1-0" (duration).
#'
#' @param n_nodes Number of nodes in the network.
#' @param directed Logical; if TRUE, the model is directed. Used for
#'   term-specific logic.
#'
#' @keywords internal
InitRedeemTerm <- function(term_name, arglist, model_type, process, n_nodes, directed = FALSE) {
  init_func_name <- paste0("InitRedeemTerm.", term_name)
  init_func <- get0(init_func_name, mode = "function", inherits = TRUE)

  if (is.null(init_func)) {
    stop(paste0("Term '", term_name, "' not recognized. No '", init_func_name, "' found."))
  }

  res <- init_func(arglist = arglist, model_type = model_type, process = process, n_nodes = n_nodes, directed = directed)


  return(res)
}

#' Check Arguments for redeem Model Terms
#'
#' Internal helper to validate arguments and set defaults.
#'
#' @param arglist List of arguments.
#' @param expected List of expected types/values.
#' @param defaults List of default values.
#' @param allowed_processes Vector of allowed processes ("0-1", "1-0").
#' @param model_type Current model type.
#'
#' @keywords internal
check.RedeemTerm <- function(arglist, expected = list(), defaults = list(),
                             allowed_processes = c("0-1", "1-0"), allowed_models = c("dem", "rem"),
                             model_type = "dem", directed = NULL,
                             directed_only = FALSE, undirected_only = FALSE) {
  # 0. Get term name for error messages
  term_name <- if (!is.null(arglist$base_name)) arglist$base_name else "unknown"

  # 1. Process check
  curr_process <- attr(arglist, "process")
  if (!is.null(curr_process) && !curr_process %in% allowed_processes) {
    stop(sprintf("Term '%s' is not allowed for the '%s' process.", term_name, curr_process))
  }

  # 2. Model check
  if (!is.null(model_type) && !model_type %in% allowed_models) {
    stop(sprintf("Term '%s' is not allowed for the '%s' model.", term_name, model_type))
  }

  # 3. Directionality check
  if (!is.null(directed)) {
    if (directed_only && !directed) {
      stop(sprintf("Term '%s' is only available for directed networks.", term_name))
    }
    if (undirected_only && directed) {
      stop(sprintf("Term '%s' is only available for undirected networks.", term_name))
    }
  }

  # 4. Defaults
  for (name in names(defaults)) {
    if (is.null(arglist[[name]])) {
      arglist[[name]] <- defaults[[name]]
    }
  }

  # 5. Expected values and types
  for (name in names(expected)) {
    val <- arglist[[name]]
    if (is.null(val)) next
    spec <- expected[[name]]
    if (is.character(spec) && length(spec) > 1) {
      if (!(val %in% spec)) {
        stop(sprintf("Argument '%s' for term '%s' must be one of: %s", name, term_name, paste(spec, collapse = ", ")))
      }
    } else if (spec == "numeric") {
      if (!is.numeric(val)) stop(sprintf("Argument '%s' for term '%s' must be numeric.", name, term_name))
    } else if (spec == "matrix") {
      if (!is.matrix(val) && !is.numeric(val)) stop(sprintf("Argument '%s' for term '%s' must be a matrix or numeric vector.", name, term_name))
    } else if (spec == "matrix_or_list") {
      if (!is.matrix(val) && !is.numeric(val) && !is.list(val)) stop(sprintf("Argument '%s' for term '%s' must be a matrix, numeric vector, or list.", name, term_name))
    } else if (spec == "numeric_or_list") {
      if (!is.numeric(val) && !is.list(val)) stop(sprintf("Argument '%s' for term '%s' must be a numeric vector or list.", name, term_name))
    }
  }

  return(arglist)
}

# --- Specific Term Initializers ---

#' @keywords internal
InitRedeemTerm.Intercept <- function(arglist, n_nodes, model_type, directed, ...) {
  check.RedeemTerm(arglist, model_type = model_type, directed = directed)
  list(base_name = "Intercept", eval_at_zero = matrix(1, n_nodes, n_nodes))
}

#' @keywords internal
InitRedeemTerm.intercept <- InitRedeemTerm.Intercept

#' @keywords internal
InitRedeemTerm.inertia <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    expected = list(
      transformation = c("identity", "log", "recip", "bin", "sig"),
      K = "numeric",
      window = "numeric"
    ),
    defaults = list(transformation = "identity", K = 1, window = Inf)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K

  base_name <- "inertia"
  list(
    base_name = base_name,
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.number_interaction <- InitRedeemTerm.inertia

#' @keywords internal
InitRedeemTerm.current_interaction <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    allowed_processes = "1-0",
    allowed_models = "dem",
    model_type = model_type,
    directed = directed,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "current_interaction", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.duration <- InitRedeemTerm.current_interaction

#' @keywords internal
InitRedeemTerm.reciprocity <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type,
    directed = directed,
    directed_only = TRUE,
    expected = list(
      transformation = c("identity", "log", "recip", "bin", "sig"),
      K = "numeric",
      window = "numeric"
    ),
    defaults = list(transformation = "identity", K = 1, window = Inf)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K

  base_name <- "reciprocity"
  list(
    base_name = base_name,
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}


#' @keywords internal
InitRedeemTerm.general_common_partner <- function(arglist, n_nodes, model_type, directed, ...) {
  if (directed) {
    arglist <- check.RedeemTerm(arglist,
      model_type = model_type, directed = directed,
      expected = list(
        transformation = c("identity", "log", "recip", "bin", "sig"),
        K = "numeric",
        type = c("OSP", "ISP", "OTP", "ITP"),
        window = "numeric"
      ),
      defaults = list(transformation = "identity", K = 1, type = "OSP", window = Inf)
    )
    base_name <- paste0("general_common_partner_", arglist$type)
  } else {
    arglist <- check.RedeemTerm(arglist,
      model_type = model_type, directed = directed,
      expected = list(
        transformation = c("identity", "log", "recip", "bin", "sig"),
        K = "numeric",
        window = "numeric"
      ),
      defaults = list(transformation = "identity", K = 1, window = Inf)
    )
    base_name <- "general_common_partner"
  }

  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(
    base_name = base_name,
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.current_common_partner <- function(arglist, n_nodes, model_type, directed, ...) {
  if (directed) {
    arglist <- check.RedeemTerm(arglist,
      model_type = model_type, directed = directed,
      allowed_models = "dem",
      expected = list(
        transformation = c("identity", "log", "recip", "bin", "sig"),
        K = "numeric",
        type = c("OSP", "ISP", "OTP", "ITP")
      ),
      defaults = list(transformation = "identity", K = 1, type = "OSP")
    )
    base_name <- paste0("current_common_partner_", arglist$type)
  } else {
    arglist <- check.RedeemTerm(arglist,
      model_type = model_type, directed = directed,
      allowed_models = "dem",
      expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
      defaults = list(transformation = "identity", K = 1)
    )
    base_name <- "current_common_partner"
  }

  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = base_name, transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.general_common_partners <- InitRedeemTerm.general_common_partner

#' @keywords internal
InitRedeemTerm.current_common_partners <- InitRedeemTerm.current_common_partner

#' @keywords internal
InitRedeemTerm.degree <- function(arglist, n_nodes, model_type, directed, ...) {
  history <- if (is.null(arglist$history)) "general" else arglist$history

  # Default type depends on whether the network is directed
  if (is.null(arglist$type)) {
    type <- if (directed) "out_sender" else "sum"
  } else {
    type <- arglist$type
  }

  is_count <- isTRUE(arglist$count)

  target_base_name <- paste0(history, "_", ifelse(is_count, "count", "degree"), "_", type)
  arglist$history <- history
  arglist$type <- type
  arglist$base_name <- target_base_name

  # Call the specific initializer
  init_func_name <- paste0("InitRedeemTerm.", target_base_name)
  init_func <- get0(init_func_name, mode = "function", inherits = TRUE)
  if (is.null(init_func)) stop(sprintf("Dispatched term '%s' not found. Check if 'history' and 'type' are valid.", target_base_name))

  init_func(arglist = arglist, n_nodes = n_nodes, model_type = model_type, directed = directed, ...)
}

#' @keywords internal
InitRedeemTerm.count <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist$count <- TRUE
  InitRedeemTerm.degree(arglist, n_nodes = n_nodes, model_type = model_type, directed = directed, ...)
}

#' @keywords internal
InitRedeemTerm.degrees <- InitRedeemTerm.degree

#' @keywords internal
InitRedeemTerm.triangle <- function(arglist, n_nodes, model_type, directed, ...) {
  history <- if (is.null(arglist$history)) "general" else arglist$history
  target_base_name <- paste0(history, "_triangle")
  arglist$history <- history
  arglist$base_name <- target_base_name

  init_func_name <- paste0("InitRedeemTerm.", target_base_name)
  init_func <- get0(init_func_name, mode = "function", inherits = TRUE)
  if (is.null(init_func)) stop(sprintf("Dispatched term '%s' not found. Check if 'history' is valid.", target_base_name))

  init_func(arglist = arglist, n_nodes = n_nodes, model_type = model_type, directed = directed, ...)
}

#' @keywords internal
InitRedeemTerm.common_partner <- function(arglist, n_nodes, model_type, directed, ...) {
  history <- if (is.null(arglist$history)) "general" else arglist$history
  target_base_name <- paste0(history, "_common_partner")
  arglist$history <- history
  arglist$base_name <- target_base_name

  init_func_name <- paste0("InitRedeemTerm.", target_base_name)
  init_func <- get0(init_func_name, mode = "function", inherits = TRUE)
  if (is.null(init_func)) stop(sprintf("Dispatched term '%s' not found. Check if 'history' is valid.", target_base_name))

  init_func(arglist = arglist, n_nodes = n_nodes, model_type = model_type, directed = directed, ...)
}

#' @keywords internal
InitRedeemTerm.common_partners <- InitRedeemTerm.common_partner

#' @keywords internal
InitRedeemTerm.general_triangle <- function(arglist, n_nodes, model_type, directed, ...) {
  if (directed) {
    arglist <- check.RedeemTerm(arglist,
      model_type = model_type, directed = directed,
      expected = list(
        transformation = c("identity", "log", "recip", "bin", "sig"),
        K = "numeric",
        type = c("OSP", "ISP", "OTP", "ITP"),
        window = "numeric"
      ),
      defaults = list(transformation = "identity", K = 1, type = "OSP", window = Inf)
    )
    base_name <- paste0("general_triangle_", arglist$type)
  } else {
    arglist <- check.RedeemTerm(arglist,
      model_type = model_type, directed = directed,
      expected = list(
        transformation = c("identity", "log", "recip", "bin", "sig"),
        K = "numeric",
        window = "numeric"
      ),
      defaults = list(transformation = "identity", K = 1, window = Inf)
    )
    base_name <- "general_triangle"
  }

  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(
    base_name = base_name,
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.current_triangle <- function(arglist, n_nodes, model_type, directed, ...) {
  if (directed) {
    arglist <- check.RedeemTerm(arglist,
      model_type = model_type, directed = directed,
      allowed_models = "dem",
      expected = list(
        transformation = c("identity", "log", "recip", "bin", "sig"),
        K = "numeric",
        type = c("OSP", "ISP", "OTP", "ITP")
      ),
      defaults = list(transformation = "identity", K = 1, type = "OSP")
    )
    base_name <- paste0("current_triangle_", arglist$type)
  } else {
    arglist <- check.RedeemTerm(arglist,
      model_type = model_type, directed = directed,
      allowed_models = "dem",
      expected = list(
        transformation = c("identity", "log", "recip", "bin", "sig"),
        K = "numeric"
      ),
      defaults = list(transformation = "identity", K = 1)
    )
    base_name <- "current_triangle"
  }

  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = base_name, transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.degree_out_sender <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "degree_out_sender", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.degree_out_receiver <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "degree_out_receiver", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.degree_in_sender <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "degree_in_sender", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.degree_in_receiver <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "degree_in_receiver", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.degree_sum <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, undirected_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "degree_sum", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.degree_absdiff <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, undirected_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "degree_absdiff", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.general_count_out_sender <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    expected = list(
      transformation = c("identity", "log", "recip", "bin", "sig"),
      K = "numeric",
      window = "numeric"
    ),
    defaults = list(transformation = "identity", K = 1, window = Inf)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K

  base_name <- "general_count_out_sender"
  list(
    base_name = base_name,
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.general_count_out_receiver <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    expected = list(
      transformation = c("identity", "log", "recip", "bin", "sig"),
      K = "numeric",
      window = "numeric"
    ),
    defaults = list(transformation = "identity", K = 1, window = Inf)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(
    base_name = "general_count_out_receiver",
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.general_count_in_sender <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    expected = list(
      transformation = c("identity", "log", "recip", "bin", "sig"),
      K = "numeric",
      window = "numeric"
    ),
    defaults = list(transformation = "identity", K = 1, window = Inf)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(
    base_name = "general_count_in_sender",
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.general_count_in_receiver <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    expected = list(
      transformation = c("identity", "log", "recip", "bin", "sig"),
      K = "numeric",
      window = "numeric"
    ),
    defaults = list(transformation = "identity", K = 1, window = Inf)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(
    base_name = "general_count_in_receiver",
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.general_count_sum <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, undirected_only = TRUE,
    expected = list(
      transformation = c("identity", "log", "recip", "bin", "sig"),
      K = "numeric",
      window = "numeric"
    ),
    defaults = list(transformation = "identity", K = 1, window = Inf)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(
    base_name = "general_count_sum",
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.general_count_absdiff <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, undirected_only = TRUE,
    expected = list(
      transformation = c("identity", "log", "recip", "bin", "sig"),
      K = "numeric",
      window = "numeric"
    ),
    defaults = list(transformation = "identity", K = 1, window = Inf)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(
    base_name = "general_count_absdiff",
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.general_degree_out_sender <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    expected = list(
      transformation = c("identity", "log", "recip", "bin", "sig"),
      K = "numeric",
      window = "numeric"
    ),
    defaults = list(transformation = "identity", K = 1, window = Inf)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K

  base_name <- "general_degree_out_sender"
  list(
    base_name = base_name,
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.general_degree_out_receiver <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    expected = list(
      transformation = c("identity", "log", "recip", "bin", "sig"),
      K = "numeric",
      window = "numeric"
    ),
    defaults = list(transformation = "identity", K = 1, window = Inf)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(
    base_name = "general_degree_out_receiver",
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.general_degree_in_sender <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    expected = list(
      transformation = c("identity", "log", "recip", "bin", "sig"),
      K = "numeric",
      window = "numeric"
    ),
    defaults = list(transformation = "identity", K = 1, window = Inf)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(
    base_name = "general_degree_in_sender",
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.general_degree_in_receiver <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    expected = list(
      transformation = c("identity", "log", "recip", "bin", "sig"),
      K = "numeric",
      window = "numeric"
    ),
    defaults = list(transformation = "identity", K = 1, window = Inf)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(
    base_name = "general_degree_in_receiver",
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.general_degree_sum <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, undirected_only = TRUE,
    expected = list(
      transformation = c("identity", "log", "recip", "bin", "sig"),
      K = "numeric",
      window = "numeric"
    ),
    defaults = list(transformation = "identity", K = 1, window = Inf)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K

  base_name <- "general_degree_sum"
  list(
    base_name = base_name,
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.general_degree_absdiff <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, undirected_only = TRUE,
    expected = list(
      transformation = c("identity", "log", "recip", "bin", "sig"),
      K = "numeric",
      window = "numeric"
    ),
    defaults = list(transformation = "identity", K = 1, window = Inf)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(
    base_name = "general_degree_absdiff",
    transformation = arglist$transformation,
    eval_at_zero = eval_at_zero,
    event_stream = arglist$event_stream,
    window = arglist$window
  )
}

#' @keywords internal
InitRedeemTerm.current_degree_out_sender <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    allowed_models = "dem", directed_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "current_degree_out_sender", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.current_degree_out_receiver <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    allowed_models = "dem", directed_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "current_degree_out_receiver", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.current_degree_in_sender <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    allowed_models = "dem", directed_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "current_degree_in_sender", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.current_degree_in_receiver <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    allowed_models = "dem", directed_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "current_degree_in_receiver", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.current_degree_sum <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    allowed_models = "dem", undirected_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "current_degree_sum", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.current_degree_absdiff <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    allowed_models = "dem", undirected_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "current_degree_absdiff", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.current_count_out_sender <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    allowed_models = "dem", directed_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "current_count_out_sender", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.current_count_out_receiver <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    allowed_models = "dem", directed_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "current_count_out_receiver", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.current_count_in_sender <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    allowed_models = "dem", directed_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "current_count_in_sender", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.current_count_in_receiver <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    allowed_models = "dem", directed_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "current_count_in_receiver", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.current_count_sum <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    allowed_models = "dem", undirected_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "current_count_sum", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.current_count_absdiff <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    allowed_models = "dem", undirected_only = TRUE,
    expected = list(transformation = c("identity", "log", "recip", "bin", "sig"), K = "numeric"),
    defaults = list(transformation = "identity", K = 1)
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- arglist$K
  list(base_name = "current_count_absdiff", transformation = arglist$transformation, eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.dyadic_cov <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    expected = list(data = "matrix_or_list", change_points = "numeric")
  )
  if (is.null(arglist$data)) stop("Term 'dyadic_cov' requires a 'data' argument.")

  data <- arglist$data
  change_points <- arglist$change_points

  if (!is.null(change_points)) {
    if (is.list(data)) {
      if (length(data) != length(change_points)) {
        stop("If 'change_points' is provided, it must have the same length as 'data'.")
      }
      names(data) <- change_points
    }
  }

  # Ensure data is a matrix of correct dimensions if it's not a scalar or a list
  if (is.list(data)) {
    if (is.null(names(data))) {
      stop("If 'data' is a list, it must be named with times or 'change_points' must be provided.")
    }
    times <- suppressWarnings(as.numeric(names(data)))
    if (any(is.na(times))) {
      stop("All names of the 'data' list must be numeric (representing times).")
    }

    # Sort by time
    ord <- order(times)
    data <- data[ord]
    times <- times[ord]

    # Validation
    for (i in seq_along(data)) {
      if (!is.matrix(data[[i]])) {
        if (length(data[[i]]) == 1) {
          data[[i]] <- matrix(data[[i]], n_nodes, n_nodes)
        } else {
          stop(sprintf("Each element of 'data' list must be a matrix or a scalar, element %d is not.", i))
        }
      }
      if (nrow(data[[i]]) != n_nodes || ncol(data[[i]]) != n_nodes) {
        stop(sprintf("Data matrix for time %s must be a %d x %d matrix.", names(data)[i], n_nodes, n_nodes))
      }
    }

    # eval_at_zero should be the snapshot with max(time <= 0)
    zero_candidates <- which(times <= 0)
    if (length(zero_candidates) == 0) {
      stop("If 'data' is a time-varying list, it must include at least one measurement time <= 0.")
    }
    zero_idx <- max(zero_candidates)
    eval_at_zero <- data[[zero_idx]]
  } else {
    if (!is.matrix(data)) {
      if (length(data) == 1) {
        data <- matrix(data, n_nodes, n_nodes)
      } else {
        stop(sprintf("Data for 'dyadic_cov' must be a matrix or a scalar, not a vector of length %d.", length(data)))
      }
    }
    if (nrow(data) != n_nodes || ncol(data) != n_nodes) {
      stop(sprintf("Data for 'dyadic_cov' must be a %d x %d matrix.", n_nodes, n_nodes))
    }
    eval_at_zero <- data
  }

  list(base_name = "dyadic_cov", transformation = "identity", data = data)
}

#' @keywords internal
InitRedeemTerm.monadic_cov <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed,
    expected = list(data = "numeric_or_list", fun = "function")
  )
  if (is.null(arglist$data)) stop("Term 'monadic_cov' requires a 'data' argument (monadic vector).")
  if (is.null(arglist$fun)) stop("Term 'monadic_cov' requires a 'fun' argument (transformation function).")

  data_monadic <- arglist$data
  fun <- arglist$fun

  # Handle time-varying monadic data (list of vectors)
  if (is.list(data_monadic)) {
    data_dyadic <- lapply(data_monadic, function(v) {
      if (length(v) != n_nodes) stop(sprintf("Monadic vector length (%d) does not match n_nodes (%d).", length(v), n_nodes))
      mat <- matrix(0, n_nodes, n_nodes)
      if (!directed) {
        # Symmetrize by applying to half and mirroring
        for (i in seq_len(n_nodes)) {
          for (j in i:n_nodes) {
            val <- fun(v[i], v[j])
            mat[i, j] <- mat[j, i] <- val
          }
        }
      } else {
        # Full matrix for directed
        for (i in seq_len(n_nodes)) {
          for (j in seq_len(n_nodes)) {
            mat[i, j] <- fun(v[i], v[j])
          }
        }
      }
      return(mat)
    })
    names(data_dyadic) <- names(data_monadic)
  } else {
    # Static monadic vector
    if (length(data_monadic) != n_nodes) stop(sprintf("Monadic vector length (%d) does not match n_nodes (%d).", length(data_monadic), n_nodes))
    data_dyadic <- matrix(0, n_nodes, n_nodes)
    if (!directed) {
      for (i in seq_len(n_nodes)) {
        for (j in i:n_nodes) {
          val <- fun(data_monadic[i], data_monadic[j])
          data_dyadic[i, j] <- data_dyadic[j, i] <- val
        }
      }
    } else {
      for (i in seq_len(n_nodes)) {
        for (j in seq_len(n_nodes)) {
          data_dyadic[i, j] <- fun(data_monadic[i], data_monadic[j])
        }
      }
    }
  }

  # Forward to dyadic_cov initializer
  arglist_dyadic <- list(data = data_dyadic, change_points = arglist$change_points)
  attr(arglist_dyadic, "process") <- attr(arglist, "process")
  res <- InitRedeemTerm.dyadic_cov(arglist_dyadic, n_nodes = n_nodes, model_type = model_type, directed = directed, ...)
  res$base_name <- "monadic_cov"
  return(res)
}

#' @keywords internal
InitRedeemTerm.baseline <- function(arglist, n_nodes, model_type, directed, ...) {
  check.RedeemTerm(arglist, model_type = model_type, directed = directed)
  # Support both singular and plural argument names for robustness, and fallback to positional 'data'
  changepoints <- if (!is.null(arglist$changepoints)) {
    arglist$changepoints
  } else if (!is.null(arglist$changepoint)) {
    arglist$changepoint
  } else {
    arglist$data
  }
  labels <- arglist$labels

  if (is.null(changepoints)) {
    stop("Term 'baseline' requires a 'changepoints' argument.")
  }

  if (!is.numeric(changepoints)) {
    stop("Argument 'changepoints' for term 'baseline' must be numeric.")
  }

  if (!is.null(labels)) {
    if (length(labels) != length(changepoints)) {
      stop("Argument 'labels' for term 'baseline' must be NULL or same length as 'changepoints'.")
    }

    # Sort together
    ord <- order(changepoints)
    changepoints <- changepoints[ord]
    labels <- labels[ord]

    # Deduplicate together (lockstep)
    dups <- duplicated(changepoints)
    changepoints <- changepoints[!dups]
    labels <- labels[!dups]
  } else {
    changepoints <- sort(unique(changepoints))
  }

  list(base_name = "baseline", changepoints = changepoints, labels = labels)
}

#' @keywords internal
InitRedeemTerm.psABBA <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    allowed_models = "rem"
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- 1
  list(base_name = "psABBA", transformation = "identity", eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.psABBY <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    allowed_models = "rem"
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- 1
  list(base_name = "psABBY", transformation = "identity", eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.psABAY <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    allowed_models = "rem"
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- 1
  list(base_name = "psABAY", transformation = "identity", eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.psABXA <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    allowed_models = "rem"
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- 1
  list(base_name = "psABXA", transformation = "identity", eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.psABXB <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    allowed_models = "rem"
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- 1
  list(base_name = "psABXB", transformation = "identity", eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.psABXY <- function(arglist, n_nodes, model_type, directed, ...) {
  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    allowed_models = "rem"
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- 1
  list(base_name = "psABXY", transformation = "identity", eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' @keywords internal
InitRedeemTerm.ps <- function(arglist, n_nodes, model_type, directed, ...) {
  mode <- if (!is.null(arglist$mode)) toupper(arglist$mode) else "ABBA"
  allowed_modes <- c("ABBA", "ABBY", "ABAY", "ABXA", "ABXB", "ABXY")
  if (!(mode %in% allowed_modes)) {
    stop(paste0("Mode '", mode, "' for term 'ps' must be one of: ", paste(allowed_modes, collapse = ", ")))
  }

  arglist <- check.RedeemTerm(arglist,
    model_type = model_type, directed = directed, directed_only = TRUE,
    allowed_models = "rem"
  )
  eval_at_zero <- matrix(0, n_nodes, n_nodes)
  eval_at_zero[1, 1] <- 1
  list(base_name = paste0("ps", mode), transformation = "identity", eval_at_zero = eval_at_zero, event_stream = arglist$event_stream)
}

#' redeem Model Terms
#'
#' @description
#' The help pages of \code{\link{rem}} and \code{\link{dem}} describe the model
#' formulation and estimation details. This page documents all statistics
#' available to be used in the model formulas, characterizing the intensities
#' of event formation and dissolution.
#'
#' In the \code{redeem} framework, models like DEM (fitted via \code{\link{dem}})
#' and REM (fitted via \code{\link{rem}}) are specified using R formulas.
#' The right-hand side of these formulas defines the structural statistics and
#' covariates, where each term must be specified separately as an explicit
#' function call (e.g., \code{~ inertia() + reciprocity(window = 10)}).
#'
#' All terms support an optional \code{transformation} argument \eqn{f}.
#' The available transformations are:
#' \itemize{
#'   \item \code{"identity"} (default): \eqn{f(x) = x}
#'   \item \code{"log"}: \eqn{f(x) = \log(x + 1)}
#'   \item \code{"recip"}: \eqn{f(x) = 1/(x+1)}
#'   \item \code{"bin"}: \eqn{f(x) = I(x > 0)}
#'   \item \code{"sig"}: sigmoid-like saturation, \eqn{f(x) = x/(x + K)}
#' }
#' Throughout, \eqn{N_{i,j}(t)} denotes the cumulative number of events from
#' \eqn{i} to \eqn{j} up to (but not including) time \eqn{t};
#' \eqn{N_{i,j}^w(t)} is the windowed analogue on \eqn{(t-w,\,t)};
#' \eqn{d_i^{\mathrm{out}}(t) = |\{l: N_{i,l}(t)>0\}|} is the historical
#' out-degree of \eqn{i}; and \eqn{c_i^{\mathrm{out}}(t) = \sum_l N_{i,l}(t)}
#' is the total event count sent by \eqn{i}.
#' The superscript \eqn{\mathrm{act}} indicates that the quantity is computed
#' on the currently active DEM network.
#'
#' The implemented terms are grouped into five categories:
#' \enumerate{
#'   \item \strong{Baseline and Nuisance Terms}: Intercept, time-varying
#'     baseline, and degree fixed effects.
#'   \item \strong{Endogenous Dyadic Terms}: Inertia, reciprocity,
#'     interaction duration, and participation shifts.
#'   \item \strong{Triadic Closure and Shared Partners}: Common partners
#'     and triangle statistics.
#'   \item \strong{Degree and Centrality Statistics}: Actor degree and
#'     event count statistics.
#'   \item \strong{Exogenous Covariates}: Dyadic and monadic covariate terms.
#' }
#'
#' @section Multi-Stream Event Covariates:
#' Most endogenous terms support covariates calculated from multiple event
#' streams. By providing an \code{event_stream} argument to a term (e.g.,
#' \code{inertia(event_stream = other_events)}), users can model one event
#' process while accounting for the history of another. The package
#' automatically handles the splintering and union of these timelines.
#'
#' @section Baseline and Nuisance Terms:
#' \itemize{
#'   \item \code{Intercept()}: Intercept: Constant log-intensity baseline.
#'     \eqn{s_{i,j}(t) = 1}. Also available as \code{intercept()}.
#'   \item \code{baseline(changepoints, labels)}: Baseline: Stepwise constant
#'     log-baseline with user-specified change points
#'     \eqn{c_1 < c_2 < \ldots < c_K}.
#'     \eqn{s_{i,j}(t) = \sum_{k=1}^{K} I(t \in [c_k, c_{k+1}))}.
#'     Coefficients are treated as nuisance parameters.
#'   \item \code{degree} / \code{degrees}: Degree Fixed Effects: Node-specific sender
#'     and receiver baselines \eqn{\alpha_i} and \eqn{\gamma_j} estimated via
#'     Minorization-Maximization (MM).
#'     Contribution to linear predictor: \eqn{\alpha_i + \gamma_j}.
#'     Treated as nuisance parameters.
#' }
#'
#' @section Endogenous Dyadic Terms:
#' \itemize{
#'   \item \code{inertia(transformation, K, event_stream, window)}: Inertia:
#'     Cumulative count of past events from \eqn{i} to \eqn{j}.
#'     \eqn{s_{i,j}(t) = f(N_{i,j}(t))}; windowed: \eqn{f(N_{i,j}^w(t))}.
#'   \item \code{reciprocity(transformation, K, event_stream, window)}: Reciprocity:
#'     Cumulative count of past events from \eqn{j} to \eqn{i}.
#'     \eqn{s_{i,j}(t) = f(N_{j,i}(t))} \strong{(Directed only)}.
#'   \item \code{current_interaction(transformation, K, event_stream)}: Duration:
#'     Time elapsed since the currently active event \eqn{(i,j)} started.
#'     \eqn{s_{i,j}(t) = f(t - t_{\mathrm{start},i,j})}
#'     \strong{(DEM only)}. Alias: \code{duration()}.
#'   \item \strong{Participation Shifts} (for two consecutive events
#'     \eqn{(A \to B) \to (C \to D)}, each statistic is 1 if the specified
#'     pattern holds, 0 otherwise; \strong{REM only, Directed only}):
#'     \itemize{
#'       \item \code{psABBA(event_stream)}: PS-ABBA: Receiver responds
#'         to sender. \eqn{s_{C,D}(t) = I(C = B,\, D = A)}.
#'       \item \code{psABBY(event_stream)}: PS-ABBY: Receiver initiates
#'         to a new target. \eqn{s_{C,D}(t) = I(C = B,\, D \ne A)}.
#'       \item \code{psABAY(event_stream)}: PS-ABAY: Sender initiates
#'         to a new target. \eqn{s_{C,D}(t) = I(C = A,\, D \ne B)}.
#'       \item \code{psABXA(event_stream)}: PS-ABXA: Outsider targets
#'         original sender. \eqn{s_{C,D}(t) = I(C \ne A, C \ne B,\, D = A)}.
#'       \item \code{psABXB(event_stream)}: PS-ABXB: Outsider targets
#'         original receiver. \eqn{s_{C,D}(t) = I(C \ne A, C \ne B,\, D = B)}.
#'       \item \code{psABXY(event_stream)}: PS-ABXY: Entirely new dyad.
#'         \eqn{s_{C,D}(t) = I(C \ne A, C \ne B,\, D \ne A, D \ne B)}.
#'       \item \code{ps(mode, event_stream)}: PS Shorthand: Dispatches to
#'         one of the six participation shift statistics above based on
#'         \code{mode} (one of \code{"ABBA"}, \code{"ABBY"}, \code{"ABAY"},
#'         \code{"ABXA"}, \code{"ABXB"}, \code{"ABXY"}).
#'     }
#' }
#'
#' @section Triadic Closure and Shared Partners:
#' \itemize{
#'   \item \code{general_common_partners(}
#'     \code{transformation, K, type,}
#'     \code{event_stream, window)}:
#'     Historical Common Partners: Number of nodes \eqn{k} sharing a historical
#'     directed path of the specified type with both \eqn{i} and \eqn{j}.
#'     \eqn{s_{i,j}(t) = f(|CP_{i,j}^{\mathrm{type}}(t)|)}.
#'     \itemize{
#'       \item \code{"OSP"} (Outgoing Shared Partner):
#'         \eqn{N_{i,k}(t)>0} and \eqn{N_{j,k}(t)>0}.
#'       \item \code{"ISP"} (Incoming Shared Partner):
#'         \eqn{N_{k,i}(t)>0} and \eqn{N_{k,j}(t)>0}.
#'       \item \code{"OTP"} (Outgoing Two-Path):
#'         \eqn{N_{i,k}(t)>0} and \eqn{N_{k,j}(t)>0}.
#'       \item \code{"ITP"} (Incoming Two-Path):
#'         \eqn{N_{k,i}(t)>0} and \eqn{N_{j,k}(t)>0}.
#'     }
#'     Aliases: \code{general_common_partner()},
#'     \code{general_common_partner_OSP()},
#'     \code{general_common_partner_ISP()},
#'     \code{general_common_partner_OTP()},
#'     \code{general_common_partner_ITP()}.
#'   \item \code{current_common_partners(}
#'     \code{transformation, K, type,}
#'     \code{event_stream)}:
#'     Active Common Partners: As \code{general_common_partners} but restricted
#'     to currently active edges. \eqn{s_{i,j}(t) =
#'     f(|CP_{i,j}^{\mathrm{type,act}}(t)|)}
#'     \strong{(DEM only)}.
#'     Aliases: \code{current_common_partner()},
#'     \code{current_common_partner_OSP()},
#'     \code{current_common_partner_ISP()},
#'     \code{current_common_partner_OTP()},
#'     \code{current_common_partner_ITP()}.
#'   \item \code{general_triangle(transformation, K, type, event_stream, window)}:
#'     Historical Triangles: Number of closed triads around \eqn{(i,j)} in the
#'     historical event network of the specified type.
#'     \eqn{s_{i,j}(t) = f(|\Delta_{i,j}^{\mathrm{type}}(t)|)}
#'     \strong{(Directed only)}.
#'   \item \code{current_triangle(transformation, K, type, event_stream)}:
#'     Active Triangles: As \code{general_triangle} but restricted to currently
#'     active edges. \strong{(DEM only, Directed only)}.
#'   \item \code{common_partner(history, type, ...)}: Common Partner Shorthand:
#'     Dispatches to \code{general_common_partners()} (\code{history="general"})
#'     or \code{current_common_partners()} (\code{history="current"}).
#'   \item \code{triangle(history, type, ...)}: Triangle Shorthand:
#'     Dispatches to \code{general_triangle()} (\code{history="general"})
#'     or \code{current_triangle()} (\code{history="current"}).
#' }
#'
#' @section Degree and Centrality Statistics:
#' \itemize{
#'   \item \code{general_degree_out_sender(}
#'     \code{transformation, K, event_stream, window)}:
#'     Sender Out-Degree: Historical out-degree of sender \eqn{i}.
#'     \eqn{s_{i,j}(t) = f(d_i^{\mathrm{out}}(t))} \strong{(Directed only)}.
#'   \item \code{general_degree_out_receiver(}
#'     \code{transformation, K, event_stream, window)}:
#'     Receiver Out-Degree: Historical out-degree of receiver \eqn{j}.
#'     \eqn{s_{i,j}(t) = f(d_j^{\mathrm{out}}(t))} \strong{(Directed only)}.
#'   \item \code{general_degree_in_sender(}
#'     \code{transformation, K, event_stream, window)}:
#'     Sender In-Degree: Historical in-degree of sender \eqn{i}.
#'     \eqn{s_{i,j}(t) = f(d_i^{\mathrm{in}}(t))} \strong{(Directed only)}.
#'   \item \code{general_degree_in_receiver(}
#'     \code{transformation, K, event_stream, window)}:
#'     Receiver In-Degree: Historical in-degree of receiver \eqn{j}.
#'     \eqn{s_{i,j}(t) = f(d_j^{\mathrm{in}}(t))} \strong{(Directed only)}.
#'   \item \code{general_degree_sum(transformation, K, event_stream, window)}:
#'     Degree Sum: Sum of historical degrees of both endpoints.
#'     \eqn{s_{i,j}(t) = f(d_i(t) + d_j(t))} \strong{(Undirected only)}.
#'   \item \code{general_degree_absdiff(}
#'     \code{transformation, K, event_stream, window)}:
#'     Degree Absolute Difference: Absolute difference in historical degrees.
#'     \eqn{s_{i,j}(t) = f(|d_i(t) - d_j(t)|)} \strong{(Undirected only)}.
#'   \item \code{general_count_out_sender(}
#'     \code{transformation, K, event_stream, window)}:
#'     Sender Out-Count: Total events sent by sender \eqn{i}.
#'     \eqn{s_{i,j}(t) = f(c_i^{\mathrm{out}}(t))} \strong{(Directed only)}.
#'   \item \code{general_count_out_receiver(}
#'     \code{transformation, K, event_stream, window)}:
#'     Receiver Out-Count: Total events sent by receiver \eqn{j}.
#'     \eqn{s_{i,j}(t) = f(c_j^{\mathrm{out}}(t))} \strong{(Directed only)}.
#'   \item \code{general_count_in_sender(}
#'     \code{transformation, K, event_stream, window)}:
#'     Sender In-Count: Total events received by sender \eqn{i}.
#'     \eqn{s_{i,j}(t) = f(c_i^{\mathrm{in}}(t))} \strong{(Directed only)}.
#'   \item \code{general_count_in_receiver(}
#'     \code{transformation, K, event_stream, window)}:
#'     Receiver In-Count: Total events received by receiver \eqn{j}.
#'     \eqn{s_{i,j}(t) = f(c_j^{\mathrm{in}}(t))} \strong{(Directed only)}.
#'   \item \code{general_count_sum(transformation, K, event_stream, window)}:
#'     Count Sum: Sum of total event counts of both endpoints.
#'     \eqn{s_{i,j}(t) = f(c_i(t) + c_j(t))} \strong{(Undirected only)}.
#'   \item \code{general_count_absdiff(}
#'     \code{transformation, K, event_stream, window)}:
#'     Count Absolute Difference: Absolute difference in total event counts.
#'     \eqn{s_{i,j}(t) = f(|c_i(t) - c_j(t)|)} \strong{(Undirected only)}.
#'   \item \code{current_degree_out_sender(transformation, K, event_stream)}:
#'     Active Sender Out-Degree: Out-degree of \eqn{i} in the active DEM network.
#'     \eqn{s_{i,j}(t) = f(d_i^{\mathrm{out,act}}(t))}
#'     \strong{(DEM only, Directed only)}.
#'   \item \code{current_degree_out_receiver(}
#'     \code{transformation, K, event_stream)}:
#'     Active Receiver Out-Degree: Out-degree of \eqn{j} in active DEM network.
#'     \eqn{s_{i,j}(t) = f(d_j^{\mathrm{out,act}}(t))}
#'     \strong{(DEM only, Directed only)}.
#'   \item \code{current_degree_in_sender(transformation, K, event_stream)}:
#'     Active Sender In-Degree: In-degree of \eqn{i} in the active DEM network.
#'     \eqn{s_{i,j}(t) = f(d_i^{\mathrm{in,act}}(t))}
#'     \strong{(DEM only, Directed only)}.
#'   \item \code{current_degree_in_receiver(}
#'     \code{transformation, K, event_stream)}:
#'     Active Receiver In-Degree: In-degree of \eqn{j} in active DEM network.
#'     \eqn{s_{i,j}(t) = f(d_j^{\mathrm{in,act}}(t))}
#'     \strong{(DEM only, Directed only)}.
#'   \item \code{current_degree_sum(transformation, K, event_stream)}:
#'     Active Degree Sum: Sum of active degrees of both endpoints.
#'     \eqn{s_{i,j}(t) = f(d_i^{\mathrm{act}}(t) + d_j^{\mathrm{act}}(t))}
#'     \strong{(DEM only, Undirected only)}.
#'   \item \code{current_degree_absdiff(transformation, K, event_stream)}:
#'     Active Degree Absolute Difference: Absolute difference in active degrees.
#'     \eqn{s_{i,j}(t) = f(|d_i^{\mathrm{act}}(t) - d_j^{\mathrm{act}}(t)|)}
#'     \strong{(DEM only, Undirected only)}.
#'   \item \code{current_count_out_sender(transformation, K, event_stream)}:
#'     Active Sender Out-Count: Total active events sent by \eqn{i}.
#'     \eqn{s_{i,j}(t) = f(c_i^{\mathrm{out,act}}(t))}
#'     \strong{(DEM only, Directed only)}.
#'   \item \code{current_count_out_receiver(}
#'     \code{transformation, K, event_stream)}:
#'     Active Receiver Out-Count: Total active events sent by \eqn{j}.
#'     \eqn{s_{i,j}(t) = f(c_j^{\mathrm{out,act}}(t))}
#'     \strong{(DEM only, Directed only)}.
#'   \item \code{current_count_in_sender(transformation, K, event_stream)}:
#'     Active Sender In-Count: Total active events received by \eqn{i}.
#'     \eqn{s_{i,j}(t) = f(c_i^{\mathrm{in,act}}(t))}
#'     \strong{(DEM only, Directed only)}.
#'   \item \code{current_count_in_receiver(}
#'     \code{transformation, K, event_stream)}:
#'     Active Receiver In-Count: Total active events received by \eqn{j}.
#'     \eqn{s_{i,j}(t) = f(c_j^{\mathrm{in,act}}(t))}
#'     \strong{(DEM only, Directed only)}.
#'   \item \code{current_count_sum(transformation, K, event_stream)}:
#'     Active Count Sum: Sum of active event counts of both endpoints.
#'     \eqn{s_{i,j}(t) = f(c_i^{\mathrm{act}}(t) + c_j^{\mathrm{act}}(t))}
#'     \strong{(DEM only, Undirected only)}.
#'   \item \code{current_count_absdiff(transformation, K, event_stream)}:
#'     Active Count Absolute Difference: Absolute difference in active event counts.
#'     \eqn{s_{i,j}(t) = f(|c_i^{\mathrm{act}}(t) - c_j^{\mathrm{act}}(t)|)}
#'     \strong{(DEM only, Undirected only)}.
#'   \item \code{degree(}
#'     \code{history, type, count, transformation, K, event_stream, window)}:
#'     Degree Shorthand: Dispatches to the appropriate degree statistic based
#'     on \code{history} (\code{"general"} or \code{"current"}) and \code{type}
#'     (\code{"out_sender"}, \code{"out_receiver"}, \code{"in_sender"},
#'     \code{"in_receiver"}, \code{"sum"}, \code{"absdiff"}).
#'     Set \code{count = TRUE} for weighted (count-based) variants.
#'     Alias: \code{degrees()}.
#'   \item \code{count(history, type, transformation, K, event_stream, window)}:
#'     Count Shorthand: Equivalent to \code{degree(..., count = TRUE)}.
#' }
#'
#' @section Exogenous Covariates:
#' \itemize{
#'   \item \code{dyadic_cov(data, change_points)}: Dyadic Covariate: Time-constant
#'     or time-varying external dyadic covariate matrix \eqn{X}.
#'     \eqn{s_{i,j}(t) = X_{i,j}(t)}.
#'   \item \code{monadic_cov(data, fun, change_points)}: Monadic Covariate: External
#'     monadic covariate vector \eqn{x} converted to a dyadic matrix via
#'     user-supplied function \eqn{g}.
#'     \eqn{s_{i,j}(t) = g(x_i(t),\, x_j(t))}.
#' }
#'
#' @param K Numeric; the evaluation point or scaling/saturation factor for the
#'   sufficient statistic (default is 1).
#' @param transformation Character; specifies the transformation to apply to the
#'   statistic. One of:
#'   \itemize{
#'     \item \code{"identity"} (default): \eqn{f(x) = x}
#'     \item \code{"log"}: \eqn{f(x) = \log(x+1)}
#'     \item \code{"recip"}: \eqn{f(x) = 1/(x+1)}
#'     \item \code{"bin"}: \eqn{f(x) = I(x>0)}
#'     \item \code{"sig"}: sigmoid-like saturation, \eqn{f(x) = x/(x+K)}
#'   }
#' @param event_stream Optional matrix or data frame; an alternative event
#'   stream to use for calculating the statistic. If \code{NULL} (default),
#'   the modeled stream is used.
#' @param window Numeric; time window for calculating the statistic (default
#'   \code{Inf}, i.e., use full history).
#' @param type Character; the specific variation of the statistic or triangle
#'   type (e.g., \code{"OSP"}, \code{"ISP"}, \code{"OTP"}, \code{"ITP"},
#'   \code{"out_sender"}, \code{"sum"}).
#' @param mode Character; the participation shift mode (e.g., \code{"ABBA"},
#'   \code{"ABBY"}).
#' @param data For \code{dyadic_cov}, a numeric matrix of dimensions
#'   \eqn{N \times N}, a scalar applied globally, or a named list of matrices
#'   for time-varying covariates. For \code{monadic_cov}, a numeric vector of
#'   length \eqn{N} or a named list of vectors for time-varying covariates.
#' @param fun A function taking two arguments \code{fun(v_i, v_j)} to generate
#'   dyadic values.
#' @param change_points Optional numeric vector; time points for time-varying
#'   covariates if \code{data} is a list.
#' @param changepoints Numeric vector; time points where the baseline intensity
#'   is allowed to change.
#' @param labels Character vector; optional labels for the resulting time
#'   intervals.
#' @param history Character; \code{"general"} for cumulative history or
#'   \code{"current"} for currently active events.
#' @param count Logical; if \code{TRUE}, uses count-based (weighted) versions
#'   of degree statistics (default \code{FALSE}).
#' @param ... Arguments passed to the underlying initialization function.
#' @return A \code{redeem_term} object (a list containing structural information about the statistic) to be used inside model formulas.
#'
#' @name redeem_terms
#' @aliases terms statistics Intercept intercept inertia reciprocity current_interaction duration number_interaction general_common_partners general_common_partner general_common_partner_OSP general_common_partner_ISP general_common_partner_OTP general_common_partner_ITP current_common_partners current_common_partner current_common_partner_OSP current_common_partner_ISP current_common_partner_OTP current_common_partner_ITP general_triangle current_triangle dyadic_cov monadic_cov baseline general_degree_out_sender general_degree_out_receiver general_degree_in_sender general_degree_in_receiver general_degree_sum general_degree_absdiff general_count_out_sender general_count_out_receiver general_count_in_sender general_count_in_receiver general_count_sum general_count_absdiff current_degree_out_sender current_degree_out_receiver current_degree_in_sender current_degree_in_receiver current_degree_sum current_degree_absdiff current_count_out_sender current_count_out_receiver current_count_in_sender current_count_in_receiver current_count_sum current_count_absdiff degree degrees count triangle common_partner psABBA psABBY psABAY psABXA psABXB psABXY ps
NULL

#' @export
#' @noRd
Intercept <- function() {
  structure(list(), class = "redeem_term", base_name = "Intercept")
}

#' @export
#' @noRd
intercept <- Intercept

#' @export
#' @noRd
inertia <- function(transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream, window = window), class = "redeem_term", base_name = "inertia")
}

#' @export
#' @noRd
number_interaction <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "number_interaction")
}

#' @export
#' @noRd
current_interaction <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "current_interaction")
}

#' @export
#' @noRd
duration <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "duration")
}

#' @export
#' @noRd
general_common_partners <- function(transformation = "identity", K = 1, type = "OSP", event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, type = type, event_stream = event_stream, window = window), class = "redeem_term", base_name = "general_common_partner")
}

#' @export
#' @noRd
general_common_partner <- general_common_partners

#' @export
#' @noRd
general_common_partner_OSP <- function(...) general_common_partners(type = "OSP", ...)
#' @export
#' @noRd
general_common_partner_ISP <- function(...) general_common_partners(type = "ISP", ...)
#' @export
#' @noRd
general_common_partner_OTP <- function(...) general_common_partners(type = "OTP", ...)
#' @export
#' @noRd
general_common_partner_ITP <- function(...) general_common_partners(type = "ITP", ...)

#' @export
#' @noRd
current_common_partners <- function(transformation = "identity", K = 1, type = "OSP", event_stream = NULL) {
  structure(list(transformation = transformation, K = K, type = type, event_stream = event_stream), class = "redeem_term", base_name = "current_common_partner")
}

#' @export
#' @noRd
current_common_partner <- current_common_partners

#' @export
#' @noRd
current_common_partner_OSP <- function(...) current_common_partners(type = "OSP", ...)
#' @export
#' @noRd
current_common_partner_ISP <- function(...) current_common_partners(type = "ISP", ...)
#' @export
#' @noRd
current_common_partner_OTP <- function(...) current_common_partners(type = "OTP", ...)
#' @export
#' @noRd
current_common_partner_ITP <- function(...) current_common_partners(type = "ITP", ...)

#' @export
#' @noRd
general_triangle <- function(transformation = "identity", K = 1, type = "OSP", event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, type = type, event_stream = event_stream, window = window),
    class = "redeem_term", base_name = "general_triangle"
  )
}

#' @export
#' @noRd
current_triangle <- function(transformation = "identity", K = 1, type = "OSP", event_stream = NULL) {
  structure(list(transformation = transformation, K = K, type = type, event_stream = event_stream),
    class = "redeem_term", base_name = "current_triangle"
  )
}

#' @export
#' @noRd
reciprocity <- function(transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream, window = window), class = "redeem_term", base_name = "reciprocity")
}

#' @export
#' @noRd
dyadic_cov <- function(data, change_points = NULL) {
  structure(list(data = data, change_points = change_points), class = "redeem_term", base_name = "dyadic_cov")
}

#' @export
#' @noRd
monadic_cov <- function(data, fun, change_points = NULL) {
  structure(list(data = data, fun = fun, change_points = change_points), class = "redeem_term", base_name = "monadic_cov")
}

#' @export
#' @noRd
baseline <- function(changepoints, labels = NULL) {
  structure(list(changepoints = changepoints, labels = labels), class = "redeem_term", base_name = "baseline")
}

#' @export
#' @noRd
general_degree_out_sender <- function(transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream, window = window), class = "redeem_term", base_name = "general_degree_out_sender")
}

#' @export
#' @noRd
general_degree_out_receiver <- function(transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream, window = window), class = "redeem_term", base_name = "general_degree_out_receiver")
}

#' @export
#' @noRd
general_degree_in_sender <- function(transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream, window = window), class = "redeem_term", base_name = "general_degree_in_sender")
}

#' @export
#' @noRd
general_degree_in_receiver <- function(transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream, window = window), class = "redeem_term", base_name = "general_degree_in_receiver")
}

#' @export
#' @noRd
general_degree_sum <- function(transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream, window = window), class = "redeem_term", base_name = "general_degree_sum")
}

#' @export
#' @noRd
general_degree_absdiff <- function(transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream, window = window), class = "redeem_term", base_name = "general_degree_absdiff")
}

#' @export
#' @noRd
general_count_out_sender <- function(transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream, window = window), class = "redeem_term", base_name = "general_count_out_sender")
}

#' @export
#' @noRd
general_count_out_receiver <- function(transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream, window = window), class = "redeem_term", base_name = "general_count_out_receiver")
}

#' @export
#' @noRd
general_count_in_sender <- function(transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream, window = window), class = "redeem_term", base_name = "general_count_in_sender")
}

#' @export
#' @noRd
general_count_in_receiver <- function(transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream, window = window), class = "redeem_term", base_name = "general_count_in_receiver")
}

#' @export
#' @noRd
general_count_sum <- function(transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream, window = window), class = "redeem_term", base_name = "general_count_sum")
}

#' @export
#' @noRd
general_count_absdiff <- function(transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream, window = window), class = "redeem_term", base_name = "general_count_absdiff")
}

#' @export
#' @noRd
current_degree_out_sender <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "current_degree_out_sender")
}

#' @export
#' @noRd
current_degree_out_receiver <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "current_degree_out_receiver")
}

#' @export
#' @noRd
current_degree_in_sender <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "current_degree_in_sender")
}

#' @export
#' @noRd
current_degree_in_receiver <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "current_degree_in_receiver")
}

#' @export
#' @noRd
current_degree_sum <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "current_degree_sum")
}

#' @export
#' @noRd
current_degree_absdiff <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "current_degree_absdiff")
}

#' @export
#' @noRd
current_count_out_sender <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "current_count_out_sender")
}

#' @export
#' @noRd
current_count_out_receiver <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "current_count_out_receiver")
}

#' @export
#' @noRd
current_count_in_sender <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "current_count_in_sender")
}

#' @export
#' @noRd
current_count_in_receiver <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "current_count_in_receiver")
}

#' @export
#' @noRd
current_count_sum <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "current_count_sum")
}

#' @export
#' @noRd
current_count_absdiff <- function(transformation = "identity", K = 1, event_stream = NULL) {
  structure(list(transformation = transformation, K = K, event_stream = event_stream), class = "redeem_term", base_name = "current_count_absdiff")
}

#' @export
#' @noRd
degree <- function(history = c("general", "current"),
                   type = c("out_sender", "out_receiver", "in_sender", "in_receiver", "sum", "absdiff"),
                   count = FALSE, transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  history <- match.arg(history)
  type <- match.arg(type)
  structure(
    list(
      history = history, type = type, count = count,
      transformation = transformation, K = K, event_stream = event_stream, window = window
    ),
    class = "redeem_term", base_name = "degree"
  )
}

#' @export
#' @noRd
degrees <- degree

#' @export
#' @noRd
count <- function(history = c("general", "current"),
                  type = c("out_sender", "out_receiver", "in_sender", "in_receiver", "sum", "absdiff"),
                  transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  history <- match.arg(history)
  type <- match.arg(type)
  degree(
    history = history, type = type, count = TRUE,
    transformation = transformation, K = K, event_stream = event_stream, window = window
  )
}

#' @export
#' @noRd
triangle <- function(history = c("general", "current"),
                     type = c("OSP", "ISP", "OTP", "ITP"),
                     transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  history <- match.arg(history)
  type <- match.arg(type)
  structure(
    list(
      history = history, type = type,
      transformation = transformation, K = K, event_stream = event_stream, window = window
    ),
    class = "redeem_term", base_name = "triangle"
  )
}

#' @export
#' @noRd
common_partner <- function(history = c("general", "current"),
                           type = c("OSP", "ISP", "OTP", "ITP"),
                           transformation = "identity", K = 1, event_stream = NULL, window = Inf) {
  history <- match.arg(history)
  type <- match.arg(type)
  structure(list(history = history, type = type, transformation = transformation, K = K, event_stream = event_stream, window = window),
    class = "redeem_term", base_name = "common_partner"
  )
}

#' @export
#' @noRd
psABBA <- function(event_stream = NULL) {
  structure(list(event_stream = event_stream), class = "redeem_term", base_name = "psABBA")
}

#' @export
#' @noRd
psABBY <- function(event_stream = NULL) {
  structure(list(event_stream = event_stream), class = "redeem_term", base_name = "psABBY")
}

#' @export
#' @noRd
psABAY <- function(event_stream = NULL) {
  structure(list(event_stream = event_stream), class = "redeem_term", base_name = "psABAY")
}

#' @export
#' @noRd
psABXA <- function(event_stream = NULL) {
  structure(list(event_stream = event_stream), class = "redeem_term", base_name = "psABXA")
}

#' @export
#' @noRd
psABXB <- function(event_stream = NULL) {
  structure(list(event_stream = event_stream), class = "redeem_term", base_name = "psABXB")
}

#' @export
#' @noRd
psABXY <- function(event_stream = NULL) {
  structure(list(event_stream = event_stream), class = "redeem_term", base_name = "psABXY")
}

#' @export
#' @noRd
ps <- function(mode = "ABBA", event_stream = NULL) {
  structure(list(mode = mode, event_stream = event_stream), class = "redeem_term", base_name = "ps")
}
