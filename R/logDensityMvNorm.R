##' Log density for the multivariate normal distrubution.
##'
##' The default for 'sigma' is the identity matrix.
##' @title Log density of the multivariate normal.
##' @param x numeric. Vector of mean values.
##' @param mu numeric. Vector of the distribution means.
##' @param sigma matrix. Variance covariance matrix parameter.
##' @return numeric. Log density.
##' @author Daniel S. Caetano, \email{caetanods1@@gmail.com}
##' @noRd
logDensityMvNorm <- function(x, mu, sigma = diag(length(x))){
    ## The log density for the multivariate normal distribution.
    ## x = vector of means.
    ## mean = true means.
    ## sigma = the variance covariance matrix. Default is an identity matrix.
    dd <- x - mu
    k <- dim(sigma)[1]
    SS <- solve(sigma)
    res <- -0.5 * (k * ( log(2) + log(pi) ) + log(det(sigma))) + (-0.5 * dd %*% SS %*% dd )
    return(res[1,1])
}
