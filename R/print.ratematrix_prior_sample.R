##' Print method for the "ratematrix_prior_sample" class.
##'
##' Print information about object.
##' @title Print method for the "ratematrix_prior_sample" class.
##' @param x The object.
##' @param ... Additional arguments. Not used here.
##' @export
print.ratematrix_prior_sample <- function(x, ...){
    ## First make some calculations:
    if( is.null( nrow( x$mu ) ) ){
        cat("\n")
        cat("1 sample from the prior distribution","\n")
        cat("\n")
        cat("Check 'names' for more details.", "\n")
    } else{
        number <- nrow( x$mu )
        cat("\n")
        cat(number, "samples from the prior distribution","\n")
        cat("\n")
        cat("Check 'names'for more details.", "\n")
    }
}
