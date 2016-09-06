##' Make the proposal and acceptance steps for the phylogenetic mean.
##'
##' Internal function to be used in the MCMC. This functions uses a simple sliding window strategy to perform the updates. This version will make independent updates for each of the phylogenetic means. This might improve the mixing.
##' @title Proposal and accept/reject for the phylogenetic mean.
##' @param cache.data list. The cache with the data.
##' @param cache.chain list. The cache with the MCMC chain.
##' @param prior list. A list with the prior functions.
##' @param v numeric. Degrees of freedom for the inverse-wishart distribution.
##' @param w_sd For standard use for the proposal of rate matrices. No use here.
##' @param w_mu Width of the proposal for the phylogenetic mean.
##' @param iter numeric. The state of the MCMC chain. This is used to access elements in the MCMC caches.
##' @param count numeric. Keep the count of the chain to record accepted and rejected steps.
##' @param traitwise Whether the proposal for the phylogenetic mean is made for all the traits or trait-by-trait.
##' @param w numeric. The sliding window width parameter.
##' @return Return a modified 'cache.chain'.
multi.phylo.mean.step.fast <- function(cache.data, cache.chain, prior, v, w_sd, w_mu, iter, count, traitwise=TRUE){

    if(traitwise == TRUE){
        select <- sample(1:cache.data$k, size=1) ## Select one of the traits to be updated.
        ## make.prop.mean is a function to make sliding window proposal moves.
        prop.root <- cache.chain$chain[[iter-1]][[1]]
        prop.root[select] <- sliding.window(prop.root[select], w_mu)
    }
    if( traitwise == FALSE){
        ## make.prop.mean is a function to make sliding window proposal moves.
        prop.root <- sapply(cache.chain$chain[[iter-1]][[1]], function(x) sliding.window(x, w_mu) )
        select <- "both (ignore NAs)"
    }

    ## Get log prior ratio. Note that the constant parameters will have a prior ratio of 1.
    prop.root.prior <- prior[[1]](prop.root)
    pp <- prop.root.prior - cache.chain$curr.root.prior

    ## Get log likelihood ratio.
    prop.root.lik <- loglikMCMC(cache.data$X, cache.data$k, cache.data$nodes, cache.data$des, cache.data$anc, cache.data$mapped.edge
                              , R=cache.chain$chain[[iter-1]][[2]], mu=as.vector(prop.root) )
    ll <-  prop.root.lik - cache.chain$lik[iter-1]
    ## Get ratio in log space.
    r <- ll + pp

    ## Acceptance step.
    ## This here need a trick on the for loop. The vcv block is the same as the nex gen.
    if(exp(r) > runif(1)){ ## Accept.
        print( "*** Root step ***" )
        print( paste0("ACCEPTED. Proposal for trait ", select, " from ", cache.chain$chain[[iter-1]][[1]][select], " to ", prop.root[select], "; r =", round(exp(r), 4), "; log_lik=", ll, "; log_prior= ", pp, ".") )
        cache.chain$chain[[iter]] <- cache.chain$chain[[iter-1]]
        cache.chain$chain[[iter]][[1]] <- prop.root
        cache.chain$curr.root.prior <- prop.root.prior
        cache.chain$acc[count] <- 1
        cache.chain$lik[iter] <- prop.root.lik
    } else{                ## Reject.
        print( "*** Root step ***" )
        print( paste0("REJECTED. Proposal for trait ", select, " from ", cache.chain$chain[[iter-1]][[1]][select], " to ", prop.root[select], "; r =", round(exp(r), 4), "; log_lik=", ll, "; log_prior= ", pp, ".") )
        cache.chain$chain[[iter]] <- cache.chain$chain[[iter-1]]
        cache.chain$acc[count] <- 0
        cache.chain$lik[iter] <- cache.chain$lik[iter-1]
    }

    ## Return cache:
    return(cache.chain)
}
