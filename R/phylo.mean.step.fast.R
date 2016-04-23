##' Make the proposal and accept/reject for the phylogenetic mean.
##'
##' Internal function to be used with the MCMC.
##' @title Phylogenetic mean proposal and accept/reject.
##' @param cache.data The cache for the data.
##' @param cache.chain The cache with the MCMC chain.
##' @param prior List with prior functions.
##' @param w The width parameter for the sliding-window proposal step of the phylogenetic mean.
##' @param v The degrees of freedom parameter for the inverted-wishart distribution.
##' @param iter Tracks the current generation of the MCMC chain.
##' @param count Used to track the accept and reject steps of the MCMC.
##' @return Updated version of the cache.chain.
phylo.mean.step.fast <- function(cache.data, cache.chain, prior, w, v, iter, count){

    ## make.prop.mean is a function to make sliding window proposal moves.
    prop.root <- sapply(cache.chain$chain[[iter-1]][[1]], function(x) sliding.window(x, w) )
    ## Get log prior ratio. Note that the constant parameters will have a prior ratio of 1.
    prop.root.prior <- prior[[1]](prop.root)
    pp <- prop.root.prior - cache.chain$curr.root.prior
    ## Create column vector format of b (phylo mean).
    #b.prop <- matrix( sapply(as.vector(prop.root), function(x) rep(x, cache.data$n) ) )
    ## Get log likelihood ratio.
    prop.root.lik <- singleR.loglik(X=cache.data$X, phy=cache.data$phy, root=as.vector(prop.root)
                                  , R=cache.chain$chain[[iter-1]][[2]], n=cache.data$n, r=cache.data$k)
    ll <-  prop.root.lik - cache.chain$lik[iter-1]
    ## Get ratio in log space.
    r <- ll + pp

    ## Acceptance step.
    ## This here need a trick on the for loop. The vcv block is the same as the nex gen.
    if(exp(r) > runif(1)){ ## Accept.
        cache.chain$chain[[iter]] <- cache.chain$chain[[iter-1]]
        #cache.chain$b.curr <- b.prop
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
