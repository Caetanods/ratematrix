##' Gerate a list of prior functions for the log density of the evolutionary rate matrix and phylogenetic mean.
##'
##' The prior for the evolutionary rate matrix uses the inverse wishart. The inverse Wishart is a very
##'     informative prior and need to be explored with caution. Please see 'make.prior.barnard' for a
##'     usually better prior option.
##' \cr
##' The prior on the phylogenetic mean (root value) is a multivariate normal distribution with mean vector
##'     'mu' and variance vector 'sd'.
##' @title Create prior distributions.
##' @param mu numeric. Mean value for the multivariate normal prior on the phylogenetic mean or root value.
##' @param sd numeric. Standard deviation for the multivariate normal prior on the phylogenetic mean.
##' @param v numeric. Degrees of freedom for the invere Wishart distribution.
##' @param p numeric. The dimension of the evolutionary rate matrix.
##' @param S matrix. The scale matrix for the inverse Wishart distribution.
##' @return List of prior functions.
##' @noRd
makePriorIWish <- function(mu, sd, v, p, S){
    ## Returns a list of the prior function for the log densitiy of the inverse-Wishart distribution
    ##      and the prior function for the phylogenetic mean vector.
    ## Prior for vcv is a inverse-Wishart centered in the given sigma with tunning parameter v.
    ## Prior for the phylogenetic mean is multivariate normal.
    ## mean = mu for the prior on the phylogenetic mean.
    ## v = the degrees of freedom of the inverse-Wishart.
    ## p = the dimension of the matrix.
    ## S = the scale matrix of the inverse-Wishart.
    prior <- list(
        prior_mu = function(x) logDensityMvNorm(x, mu, sigma = diag(sd^2, length(x))),
        #prior_S = function(x) logDensityIWish(x, v = v, S = ((v-p-1) * S) )
		## Why we have '(v-p-1)' in this prior function. This seems strange.
		prior_S = function(x) logDensityIWish(x, v = v, S = ((v-p+1) * S) )
		## Note the change above from the '-' to the '+'. This is following Barnard.
        )
    return(prior)
}
