##' Makes the proposal and acceptance step for the R matrix.
##'
##' Internal function to be used with the MCMC.
##' @title Proposal and accept reject of the evolutionary rate matrix.
##' @param cache.data Cache of data.
##' @param cache.chain Cache of the MCMC chain
##' @param prior List with the prior functions
##' @param w Width of the sliding-window parameter
##' @param v Degrees of freedom parameter for the inverted-wishart.
##' @param iter Records the generation of the MCMC chain.
##' @param count Keep track of the accept reject.
##' @return Updated cache.chain.
R.matrix.step.fast <- function(cache.data, cache.chain, prior, w, v, iter, count){
    ## Makes the proposal and acceptance step for the R matrix.
    ## cache.data = cache with data for analysis.
    ## cache.chain = cache with chain objects.
    ## prior = prior functions.
    ## w = sliding window width parameter.
    ## v = degree of freedom parameter of the inverse-Wishart.
    ## iter = the state of the MCMC chain. This is used to access data in the cache.chain
    ##        objects.

    ## make.prop.vcv is a function to make a Barnard separation proposal for vcv.
    prop.vcv <- make.prop.iwish(cache.chain$chain[[iter-1]][[2]], k=cache.data$k, v=v)
    ## Calculates the hastings for the inverse-Wishart. Uses the ration of transition
    ##    densities only.
    hh <- h.density(cache.chain$chain[[iter-1]][[2]], prop.vcv, p=cache.data$k, v=v)
    ## Get log prior ratio. Note that the constant parameters will have a prior ratio of 1.
    prop.vcv.prior <- prior[[2]](prop.vcv)
    pp <- prop.vcv.prior - cache.chain$curr.vcv.prior
    ## Get log likelihood ratio.
    prop.vcv.lik <- singleR.loglik(X=cache.data$X, root=as.vector(cache.chain$chain[[iter-1]][[1]])
                                       , R=prop.vcv, C=cache.data$C, D=cache.data$D
                                       , n=cache.data$n, r=cache.data$k, method="rpf")
    ## prop.vcv.lik <- singleR.loglik(R=prop.vcv, C=cache.data$C, y=cache.data$y, b=cache.chain$b.curr
    ##                              , n=cache.data$n, r=cache.data$k)
    ll <-  prop.vcv.lik - cache.chain$lik[iter-1]
    ## Get ratio in log space.
    r <- ll + pp + hh

    ## Acceptance step.
    ## This here need a trick on the for loop. The vcv block is the same as the nex gen.
    if(exp(r) > runif(1)){ ## Accept.
        cache.chain$chain[[iter]] <- cache.chain$chain[[iter-1]]
        cache.chain$chain[[iter]][[2]] <- prop.vcv
        cache.chain$curr.vcv.prior <- prop.vcv.prior
        cache.chain$acc[count] <- 2
        cache.chain$lik[iter] <- prop.vcv.lik
    } else{                ## Reject.
        cache.chain$chain[[iter]] <- cache.chain$chain[[iter-1]]
        cache.chain$acc[count] <- 0
        cache.chain$lik[iter] <- cache.chain$lik[iter-1]
    }

    ## Return cache:
    return(cache.chain)
}
