##' make a proposal for a new covariance matrix given the current state.
##'
##' Internal function to be used in the MCMC.
##' @title Proposal for covariance matrices.
##' @param curr.vcv matrix. Current state of the covariance matrix.
##' @param k numeric. Number of dimensions of the matrix.
##' @param v numeric. Degrees of freedom for the inverse wishart distribution.
##' @return A matrix.
##' @noRd
makePropIWish <- function(curr.vcv, k, v){
    ## Make a proposal for a new vcv matrix given the current state.
    ## k = dimension of the matrix.
    ## v = degrees of freedom.
    center <- (v-k-1) * curr.vcv
    prop <- riwish(v, center)
    return(prop)
}
