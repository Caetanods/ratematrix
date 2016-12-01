##' Gerate a list of prior functions for the log density of the evolutionary rate matrix and phylogenetic mean.
##'
##' The prior for the evolutionary rate matrix uses the separation strategy (Barnard et al, 2000).
##'     This strategy creates a marginally uniform prior by modelling the covariance matrix using a
##'     correlation matrix and a vector of variance coeficients. One can also use a inverse Wishart prior,
##'     however the inverse Wishart is a very informative prior and need to be explored with caution.\cr
##' \cr
##' The prior on the phylogenetic mean (root value) is a multivariate normal distribution with mean vector
##'     'mu' and variance vector 'sd'.
##' @title Create prior distributions.
##' @param mu numeric. Mean value for the multivariate normal prior on the phylogenetic mean or root value.
##' @param sd numeric. Standard deviation for the multivariate normal prior on the phylogenetic mean.
##' @param min numeric. Minimum bound of the uniform prior on the variance component of the covariance matrix.
##' @param max numeric. Maximum bound of the uniform prior on the variance component of the covariance matrix.
##' @return List of prior functions..
##' @references Barnard, J., R. McCulloch, and X.-L. Meng. 2000. Modeling covariance matrices in terms of standard deviations and correlations, with application to shrinkage. Statistica Sinica 10:1281–1312.
##' @export
makePriorBarnard <- function(mu, sd, min, max){
    ## Returns a list of the prior function for the log density of the vcv matrix using the Barnard
    ##      separation strategy and the multivariate normal for the root vector.
    ## Prior for the vcv is a Barnard strategy. With an invert-Wishart for the correlation matrix
    ##      and a uniform for the variance component.
    ## Prior for the phylogenetic mean is a multivariate normal.
    ## mu = mean for the phylogenetic mean.
    ## sd = standard deviation for the phylogenetic mean.
    ## min and max = bounds for the uniform distribution on the variance component of the
    ##      vcv matrix.
    prior <- list(
        prior_mu = function(x) logDensityMvNorm(x, mu, sigma = diag(sd^2, length(x))),
	prior_S = function(x) priorDensityVcvBarnard(x, min=min, max=max)
        )
    return(prior)
}
