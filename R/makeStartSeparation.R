##' Create a list with the starting points for the MCMC.
##'
##' Create a list with the starting points for the MCMC.
##' @title Start point for MCMC.
##' @param n.traits Number of traits in the model.
##' @param n.R Number of R matrices fitted to the tree.
##' @return A list with the starting points for a MCMC. This function will draw from an uniform prior. Same as created by 'make.prior.zhang'.
##' @export
##' @importFrom corpcor decompose.cov
makeStartSeparation <- function(n.traits, n.R){
    mn <- runif(n.traits, min=-100, max=100)
    if( n.R == 1 ){
        corr <- riwish(v=n.traits+1, S=diag(nrow=n.traits))
        sd <- sqrt( decompose.cov(corr)$v )
    } else{
        corr <- list()
        sd <- list()
        for(i in 1:n.R){
            corr[[i]] <- riwish(v=n.traits+1, S=diag(nrow=n.traits))
            sd[[i]] <- sqrt( decompose.cov(corr[[i]])$v )
        }
    }
    start <- list(mn, corr, sd)
    return(start)
}
