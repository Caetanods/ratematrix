##' Make the proposal and accept/reject for the phylogenetic mean.
##'
##' Internal function to be used with the MCMC.
##' @title Phylogenetic mean proposal and accept/reject.
##' @param cache.data The cache for the data.
##' @param cache.chain The cache with the MCMC chain.
##' @param prior List with prior functions.
##' @param w_sd 
##' @param w_mu 
##' @param v The degrees of freedom parameter for the inverted-wishart distribution.
##' @param iter Tracks the current generation of the MCMC chain.
##' @param count Used to track the accept and reject steps of the MCMC.
##' @param traitwise 
##' @param w The width parameter for the sliding-window proposal step of the phylogenetic mean.
##' @return Updated version of the cache.chain.
makePropMean <- function(cache.data, cache.chain, prior, w_sd, w_mu, v, iter, count, traitwise){

    if(traitwise == TRUE){
        select <- sample(1:cache.data$k, size=1) ## Select one of the traits to be updated.
        ## make.prop.mean is a function to make sliding window proposal moves.
        prop.root <- cache.chain$chain[[iter-1]][[1]]
        prop.root[select] <- slideWindow(prop.root[select], w_mu)
    }
    if(traitwise == FALSE){
        ## make.prop.mean is a function to make sliding window proposal moves.
        prop.root <- sapply(cache.chain$chain[[iter-1]][[1]], function(x) slideWindow(x, w_mu) )
        ## select <- "both (ignore NAs)" ## This is to write the accept reject to the log. Not implemented for the single matrix case.
    }

    ## make.prop.mean is a function to make sliding window proposal moves.
    ## prop.root <- sapply(cache.chain$chain[[iter-1]][[1]], function(x) slideWindow(x, w_mu) )
    ## Get log prior ratio. Note that the constant parameters will have a prior ratio of 1.
    prop.root.prior <- prior[[1]](prop.root)
    pp <- prop.root.prior - cache.chain$curr.root.prior
    ## Create column vector format of b (phylo mean).
    #b.prop <- matrix( sapply(as.vector(prop.root), function(x) rep(x, cache.data$n) ) )
    ## Get log likelihood ratio.
    prop.root.lik <- logLikSingleRegime(data=cache.data, chain=cache.chain, root=as.vector(prop.root)
                                  , R=cache.chain$chain[[iter-1]][[4]] )
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
