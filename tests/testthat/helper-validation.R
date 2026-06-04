#' Helper to close any open interactions in an event matrix
#' 
#' This is useful for simulation tests where the process stops abruptly,
#' leaving some interactions open. Mandatory validation in dem() requires
#' all interactions to be closed.
