##' Log density of the scaled inverse chi-squared distribution
##'
##' Calculates the density of the scaled inverse chi-squared distribution in log space. \cr
##' \cr
##' The parameter 'scale' has a default value of 1/v.chi.
##' @title Log density of the scaled inverse chi-squared
##' @param x vector. Vector of variances.
##' @param v.chi numeric. Degrees of freedom.
##' @param scale numeric. Scale for the distribution.
##' @return numeric. Log density of the scaled inverse chi-squared.
##' @author Daniel S. Caetano, \email{caetanods1@@gmail.com}
dinvchisq <- function(x, v.chi, scale = 1/v.chi){
    ## Returns the log density of the scaled inverse chi-squared distribution.
    ## x = vector of variances.
    ## v.chi = degrees of freedom.
    nu <- v.chi/2
    d <- nu * log(nu) - log(gamma(nu)) + nu * log(scale) - (nu + 1) * log(x) - (nu * scale/x)
    return(d)
}
