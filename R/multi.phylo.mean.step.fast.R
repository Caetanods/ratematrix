##' Make the proposal and acceptance steps for the phylogenetic mean.
##'
##' Internal function to be used in the MCMC. This functions uses a simple sliding window strategy to perform the updates.
##' @title Proposal and accept/reject for the phylogenetic mean.
##' @param cache.data list. The cache with the data.
##' @param cache.chain list. The cache with the MCMC chain.
##' @param prior list. A list with the prior functions.
##' @param w numeric. The sliding window width parameter.
##' @param v numeric. Degrees of freedom for the inverse-wishart distribution.
##' @param iter numeric. The state of the MCMC chain. This is used to access elements in the MCMC caches.
##' @param count numeric. Keep the count of the chain to record accepted and rejected steps.
##' @return Return a modified 'cache.chain'.
multi.phylo.mean.step.fast <- function(cache.data, cache.chain, prior, v, w_sd, w_mu, iter, count){
    ## Make the proposal and accept step for the phylogenetic mean.
    ## cache.data = cache with data for analysis.
    ## cache.chain = cache with chain objects.
    ## prior = prior functions.
    ## w = sliding window width parameter.
    ## v = degree of freedom parameter of the inverse-Wishart.
    ## iter = the state of the MCMC chain. This is used to access data in the cache.chain
    ##        objects.

    ## make.prop.mean is a function to make sliding window proposal moves.
    prop.root <- sapply(cache.chain$chain[[iter-1]][[1]], function(x) sliding.window(x, w_mu) )
    ## Get log prior ratio. Note that the constant parameters will have a prior ratio of 1.
    prop.root.prior <- prior[[1]](prop.root)
    pp <- prop.root.prior - cache.chain$curr.root.prior
    ## Create column vector format of b (phylo mean).
    ##b.prop <- matrix( sapply(as.vector(prop.root), function(x) rep(x, cache.data$n) ) )
    ## Get log likelihood ratio.
    prop.root.lik <- loglikMCMC(cache.data$X, cache.data$k, cache.data$nodes, cache.data$des, cache.data$anc, cache.data$mapped.edge
                              , R=cache.chain$chain[[iter-1]][[2]], mu=as.vector(prop.root) )
    ll <-  prop.root.lik - cache.chain$lik[iter-1]
    ## Get ratio in log space.
    r <- ll + pp

    ## Acceptance step.
    ## This here need a trick on the for loop. The vcv block is the same as the nex gen.
    if(exp(r) > runif(1)){ ## Accept.
        cache.chain$chain[[iter]] <- cache.chain$chain[[iter-1]]
        ## cache.chain$b.curr <- b.prop
        cache.chain$chain[[iter]][[1]] <- prop.root
        cache.chain$curr.root.prior <- prop.root.prior
        cache.chain$acc[count] <- 1
        cache.chain$lik[iter] <- prop.root.lik
    } else{                ## Reject.
        cache.chain$chain[[iter]] <- cache.chain$chain[[iter-1]]
        cache.chain$acc[count] <- 0
        cache.chain$lik[iter] <- cache.chain$lik[iter-1]
    }

    ## Return cache:
    return(cache.chain)
}
