##' Print method for the "ratematrix_prior_sample" class.
##'
##' Print information about object.
##' @title Print method for the "ratematrix_prior_sample" class.
##' @param x The object.
##' @export
print.ratematrix_prior_sample <- function(x){
    ## First make some calculations:
    number <- nrow( x$mu )
    cat(number, "samples from the prior distribution","\n")
    cat("\n")
    cat("Check 'names' of x for more details.", "\n")
}
