##' Print method for the "ratematrix_prior_function" class.
##'
##' Print information about object.
##' @title Print method for the "ratematrix_prior_function" class.
##' @param x The object.
##' @param ... Additional arguments. Not used here.
##' @export
print.ratematrix_prior_function <- function(x, ...){
    cat("\n")
    cat("Prior function for the MCMC with", as.numeric(x$pars$r), "traits and", as.numeric(x$pars$p),"regimes.", "\n")
}
