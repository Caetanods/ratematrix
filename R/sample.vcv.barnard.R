##' Sample a variance covariance matrix from the prior.
##'
##' This function samples a variance covariance matrix from the prior distribution using the Barnard et al (2009) strategy to separate the matrix in its variance and correlation component. This function is used to initiate a Markov chain Monte Carlo run by sampling a evolutionary rate matrix from its prior. 
##' @title Samples from the covariance prior.
##' @param k numeric. dimension of the matrix.
##' @param min numeric. minimum value for the uniform prior on the variance.
##' @param max numeric. maximum value for the uniform prior on the variance.
##' @return Returns a covariance matrix.
##' @export
sample.vcv.barnard <- function(k, min=0, max=100){
    ## Function to sample from the prior. Using the separation strategy of Barnard.
    ## k = the dimension of the matrix.
    ## min = minimum for the uniform prior on the variance.
    ## max = maximum for the uniform prior on the variance.
    ## Note that 'min' and 'max' need to have only positive values.
    
    ## To get the correlation matrix we need to sample from an inverse-wishart and
    ##      transform the result in a correlation matrix. (Scale matrix is I)
    S <- riwish(v=k+1, S=diag(k))
    R <- cov2cor(S)
    #var <- rinvchisq(n = k, v.chi = v.chi, scale = 1/v.chi) #Not used. The inverse scaled chisquared.
    var <- runif(n = k, min=0, max=100)
    vcv <- rebuild.cov(r=R, v=var)
    return(vcv)
}
