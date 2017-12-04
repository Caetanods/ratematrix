##' Make the proposal and accept/reject for the phylogenetic mean.
##'
##' Internal function to be used with the MCMC.
##' @title Phylogenetic mean proposal and accept/reject.
##' @param cache.data The cache for the data.
##' @param cache.chain The cache with the MCMC chain.
##' @param prior List with prior functions.
##' @param w_sd the width for the vector of standard deviations
##' @param w_mu the width for the root values
##' @param v The degrees of freedom parameter for the inverted-wishart distribution.
##' @param iter Tracks the current generation of the MCMC chain.
##' @param count Used to track the accept and reject steps of the MCMC.
##' @param w The width parameter for the sliding-window proposal step of the phylogenetic mean.
##' @return Updated version of the cache.chain.
##' @noRd
makePropMean <- function(cache.data, cache.chain, prior, w_sd, w_mu, v, files, phy){

    ## make.prop.mean is a function to make sliding window proposal moves.
    to.update.mu <- cache.chain$chain[[1]]
    prop.root <- sapply(1:length(to.update.mu), function(x) slideWindow(to.update.mu[x], w_mu[x]) )
    ## select <- "both (ignore NAs)" ## This is to write the accept reject to the log. Not implemented for the single matrix case.

    ## make.prop.mean is a function to make sliding window proposal moves.
    ## prop.root <- sapply(cache.chain$chain[[iter-1]][[1]], function(x) slideWindow(x, w_mu) )
    ## Get log prior ratio. Note that the constant parameters will have a prior ratio of 1.
    prop.root.prior <- prior[[1]](prop.root)
    pp <- prop.root.prior - cache.chain$curr.root.prior
    ## Create column vector format of b (phylo mean).
    #b.prop <- matrix( sapply(as.vector(prop.root), function(x) rep(x, cache.data$n) ) )
    ## Get log likelihood ratio.
    prop.root.lik <- logLikSingleRegime(data=cache.data, chain=cache.chain, phy=phy, root=as.vector(prop.root)
                                      , R=cache.chain$chain[[4]] )
    ll <-  prop.root.lik - cache.chain$lik
    ## Get ratio in log space.
    r <- ll + pp

    ## Acceptance step.
    ## This here need a trick on the for loop. The vcv block is the same as the nex gen.
    if(exp(r) > stats::runif(1)){ ## Accept.
        cat("1; 0; 0; 1; 1; ", prop.root.lik, "\n", sep="", file=files[[2]], append=TRUE)
        cache.chain$chain[[1]] <- prop.root
        cache.chain$curr.root.prior <- prop.root.prior
        cache.chain$lik <- prop.root.lik
    } else{                ## Reject.
        cat("0; 0; 0; 1; 1; ", cache.chain$lik, "\n", sep="", file=files[[2]], append=TRUE)
    }

    ## Return cache:
    return(cache.chain)
}
