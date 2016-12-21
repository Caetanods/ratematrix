##' Calculates the log density of the inverse Wishart distribution.
##'
##' Modified from the 'diwish' function of the package 'MCMCpack'.
##' @title Log density for the inverse Wishart.
##' @param W matrix. The covariance matrix.
##' @param v numeric. The degrees of freedom parameter.
##' @param S matrix. The standard matrix for the distribution.
##' @return numeric. The log density.
##' @author Daniel S. Caetano, \email{caetanods1@@gmail.com}
##' @noRd
logDensityIWish <- function(W, v, S){
    ## This function is derived from the 'diwish' from MCMCpack.
    ## This returns the log density.
    k <- nrow(S)
    lgammapart <- 0
    for (i in 1:k) {
        lgammapart <- lgammapart + lgamma((v + 1 - i)/2)
    }
    ldenom <- lgammapart + ( ( (v*k)/2 ) * log( 2 ) ) + ( ( (k*(k-1))/4 ) * log( pi ) )
    detS <- det(S)
    detW <- det(W)
    hold <- S %*% solve(W)
    tracehold <- sum(hold[row(hold) == col(hold)])
    lnum <- ( ( v/2 ) * log( detS ) ) + ( ( -(v + k + 1)/2 ) * log( detW ) ) + ( -1/2 * tracehold )
    return(lnum - ldenom)
}
