##' A simple sliding window proposal.
##'
##' A simple sliding window proposal. Internal function to be used by the MCMC.
##' @title Sliding window proposal.
##' @param x numeric. The current value.
##' @param w numeric. The width of the proposal.
##' @return A updated value based on 'x'.
##' @noRd
slideWindow <- function(x, w){
    ## Sliding window proposal for unbounded trait.
    ## x = the current value.
    ## w = the width parameter of the proposal.
    y <- stats::runif(1, min = x - (w/2), max = x + (w/2) )
    return(y)
}
