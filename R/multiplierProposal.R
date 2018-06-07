##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Multiplier proposal
##' @param x 
##' @param a 
##' @return The proposal and the proposal ratio.
##' @author daniel
##' @noRd
##' @importFrom stats runif setNames
multiplierProposal <- function(x, a){
    ## This proposal scheme will perform steps in the log space of the parameter.
    ## This is much more efficient for the cases in which small changes in large
    ##    values effect less than small changes in small values.
    ## Note this proposal will never change the sign of the parameter.
    ## x = the current value.
    ## a = the scale of the multiplier.
    lambda <- 2 * log(a)
    m <- exp( lambda * runif(1, min=-0.5, max=0.5) )
    ## Note that 'm' here is the proposal ratio. So need to spit this out.
    return( setNames(c(m * x, m), c("prop","prop.ratio") ) )
}
