##' Function to calculate the prior density using the Barnard separation.
##'
##' Internal function to be used by the MCMC.
##' @title Prior density under the Barnard separation.
##' @param S The variance covariance matrix. The evolutionary rate matrix.
##' @param min Minimum of the uniform prior on the variance component.
##' @param max Maximum of the uniform prior on the variance component.
##' @return The log density of the matrix S under the prior.
##' @importFrom corpcor decompose.cov
##' @noRd
priorDensityVcvBarnard <- function(S, min=0, max=100){
    ## Function to calculate the prior density using the Barnard separation.
    ## Here both the variance and the correlations are uniform. The iwish for the
    ##      correlation matrix is set to k+1, which yelds a uniform distribution for
    ##      the correlations.
    ## The variance component has a uniform distribution.
    ## S = a variance-covariance matrix.
    ## min = the minimum for the uniform prior on the variance.
    ## max = the maximum for the uniform prior on the variance.
    k <- dim(S)[1]
    ## 'decompose.cov' gets the correlation and the variance from a vcv.
    ## 'rebuild.cov' reconstruct the vcv from the correlation and the variance vector.
    decom <- decompose.cov(S)
    R <- decom$r
    var <- decom$v
    if(k == 2){
        sub.prod <- prod( t( sapply(1:k, function(x) R[-x,-x] ) ) )
        ## Here the det() of a single value is equal to the value.
    } else{
        sub.prod <- prod( t( sapply(1:k, function(x) det(R[-x,-x]) ) ) )
    }
    dR <- (log(det(R))*(k*(k-1)/2-1)) + (log(sub.prod)*(-(k+1)/2))
    ## There is a more general expression in Barnard. If informative priors are needed.
    dvar <- sum( stats::dunif(x = var, min=min, max=max, log=TRUE) ) ## Indepedent and Identically Distributed.
    d <- dR + dvar
    ## This follows Barnand and assumes the variance component is independent of the
    ##      correlation.
    return(d)
}
