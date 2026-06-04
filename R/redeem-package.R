#' redeem: Scalable Durational Event Models
#'
#' The \pkg{redeem} package provides a comprehensive framework for the estimation
#' of Durational Event Models (DEM) and Relational Event Models (REM).
#'
#' @details
#' The package is built to handle relational event sequences where interactions
#' have specific start and end times. It implements the scalable block-coordinate
#' ascent algorithm described in Fritz et al. (2026), which allows for the
#' estimation of high-dimensional actor-specific effects and time-varying baseline
#' intensities.
#'
#' Key functions:
#' \itemize{
#'   \item \code{\link{dem}}: Estimates a Durational Event Model.
#'   \item \code{\link{rem}}: Estimates a Relational Event Model.
#' }
#'
#' @author Cornelius Fritz
#'
#' @references
#' Fritz, C., Rastelli, R., Fop, M., & Caimo, A. (2026). Scalable Durational Event Models:
#' Application to Physical and Digital Interactions. arXiv:2504.00049.
#'
#' @importFrom stats binomial poisson
#' @name redeem
#' @keywords internal
"_PACKAGE"

# Suppress R CMD check NOTEs for data.table column references used in
# non-standard evaluation (`:=` assignments and `[i, j]` expressions).
utils::globalVariables(c("offset", "diff", "event", "from_avail", "to_avail", "avail"))
