##' Generate samples from an scaled inverse chi-squared distribution.
##'
##' Internal function to be used by the MCMC.
##' @title Random sample from scaled inverse chi-squared distribution.
##' @param n numeric. Number of samples.
##' @param v.chi numeric. Degrees of freedom.
##' @param scale numeric. The scale factor.
##' @return numeric. Sample from the distribution.
rinvchisq <- function(n = 1, v.chi = 1, scale = 1/v.chi){
    ## Sample from the scaled inverse chi-squared distribution.
    ## n = number of samples.
    ## v.chi = degrees of freedom.
    return( (v.chi * scale) / rchisq(n, df = v.chi) )
}
