.onAttach <- function(lib, pkg) {
	packageStartupMessage("Welcome to GWmodel version 2.0-8.\n Note: This verision has been re-built with RcppArmadillo to improve its performance.\nThe new version of GWmodel 2.0-8 now is ready, and new functions have been incorporated or improved:\n 1) gwr.multiscale, multiscale GWR\n", appendLF = FALSE)
  }