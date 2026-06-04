#' @useDynLib redeem, .registration=TRUE

.onUnload <- function (libpath) {
  library.dynam.unload("redeem", libpath)
}
