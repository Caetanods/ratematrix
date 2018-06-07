##' Print method for the "ratematrix_multi_chain" class.
##'
##' Print information about object.
##' @title Print method for the "ratematrix_multi_chain" class.
##' @param x The object.
##' @param ... Additional arguments. Not used here.
##' @export
print.ratematrix_multi_chain <- function(x, ...){
    ## First make some calculations:
    k <- ncol( x$root )
    post.length <- nrow( x$root )
    p <- length( x$matrix )
    cat("\n")
    cat("Posterior distribution with multiple regimes","\n")
    cat("Number of traits: ", k, "\n")
    cat("Number of regimes: ", p, "\n")
    cat("Number of posterior samples: ", post.length, "\n")
    cat("\n")
    cat("Use 'plotRatematrix' and 'plotRootValue' to plot the distribution.", "\n")
    cat("Use 'checkConvergence' to verify convergence.", "\n")
    cat("Use 'mergePosterior' to merge two or more posterior chains.", "\n")
    cat("Check 'names' for more details.", "\n")
}
