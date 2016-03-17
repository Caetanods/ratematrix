##' Makes the proposal and accept/reject steps for two or more fitted evolutionary rate matrices.
##'
##' Internal function to be used by the MCMC
##' @title Proposal and accept/reject for two or more evolutionary rate matrices.
##' @param cache.data list. Cache with the data for the MCMC
##' @param cache.chain list. Cache with the states of the MCMC chain
##' @param prior list. List of prior functions.
##' @param w numeric. Width of the sliding-window proposal for the phylogenetic mean.
##' @param v numeric. Degrees of freedom for the inverse-Wishart distribution.
##' @param iter numeric. Gives the present state (generation) of the MCMC chain.
##' @param count numeric. Used to track the accept and reject steps of the MCMC chain.
##' @return Updated state of the cache.chain list.
multiR.matrix.step.fast <- function(cache.data, cache.chain, prior, w, v, iter, count){
    ## Makes the proposal and acceptance step for the multiR matrix MCMC.
    ## This version assumes the prior function is the same for all R matrices.
    ## Also applies the same proposal step.
    ## Beware that some objects are now assumed to be lists.
    ## cache.data = cache with data for analysis.
    ## cache.chain = cache with chain objects.
    ## prior = prior functions.
    ## w = sliding window width parameter.
    ## v = degree of freedom parameter of the inverse-Wishart.
    ## iter = the state of the MCMC chain. This is used to access data in the cache.chain
    ##        objects.

    ## Random draw one of the p R matrices to update:
    Rp <- sample(1:cache.data$p, size=1)

    ## make.prop.vcv is a function to make a Barnard separation proposal for vcv.
    ## Note that 'Rp' points to one of the p R matrices.
    prop.vcv <- make.prop.iwish(cache.chain$chain[[iter-1]][[2]][[Rp]], k=cache.data$k, v=v)
    ## Calculates the hastings for the inverse-Wishart. Uses the ratio of transition
    ##    densities only.
    ## The hastings is only based on the proposed R matrix, since the rest is constant and
    ##    have hastings ratio equal to 1.
    hh <- h.density(cache.chain$chain[[iter-1]][[2]][[Rp]], prop.vcv, p=cache.data$k, v=v)
    ## Get log prior ratio. Note that the constant parameters will have a prior ratio of 1.
    prop.vcv.prior <- prior[[2]](prop.vcv)
    pp <- prop.vcv.prior - cache.chain$curr.vcv.prior[[Rp]]
    ## Get log likelihood ratio.
    ## Need to create 'prop.vcv' object with all p R matrices.
    prop.vcv.list <- cache.chain$chain[[iter-1]][[2]] ## List with size p.
    prop.vcv.list[[Rp]] <- prop.vcv ## Change only the proposed R matrix.
    #prop.vcv.lik <- multiR.loglik(R.m=prop.vcv.list, C.m=cache.data$C.m, y=cache.data$y, b=cache.chain$b.curr
                                  #, n=cache.data$n, r=cache.data$k, p=cache.data$p)
    prop.vcv.lik <- multiR.loglik(X=cache.data$X, root=as.vector(cache.chain$chain[[iter-1]][[1]])
                                , R.m=prop.vcv.list, C.m=cache.data$C.m, D=cache.data$D
                                , n=cache.data$n, r=cache.data$k, p=cache.data$p, method="rpf")
    ll <-  prop.vcv.lik - cache.chain$lik[iter-1]
    ## Get ratio in log space.
    r <- ll + pp + hh

    ## Acceptance step.
    ## This here need a trick on the for loop. The vcv block is the same as the nex gen.
    
    if(exp(r) > runif(1)){ ## Accept.
        cache.chain$chain[[iter]] <- cache.chain$chain[[iter-1]]
        ## Update only one of the R matrices.
        cache.chain$chain[[iter]][[2]][[Rp]] <- prop.vcv
        cache.chain$curr.vcv.prior[[Rp]] <- prop.vcv.prior
        ## The index for the acceptance is 'Rp+1', so we can track which of the matrices are accepted.
        cache.chain$acc[count] <- Rp+1
        cache.chain$lik[iter] <- prop.vcv.lik
    } else{                ## Reject.
        cache.chain$chain[[iter]] <- cache.chain$chain[[iter-1]]
        cache.chain$acc[count] <- 0
        cache.chain$lik[iter] <- cache.chain$lik[iter-1]
    }

    ## Return cache:
    return(cache.chain)
}
