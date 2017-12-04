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
##' @param files the path to write to the log file.
##' @return Return a modified 'cache.chain'.
##' @importFrom MASS mvrnorm
##' @noRd
makePropMeanForMult <- function(cache.data, cache.chain, prior, v, w_sd, w_mu, files){

    ## make.prop.mean is a function to make sliding window proposal moves.
    to.update.mu <- cache.chain$chain[[1]]
    prop.root <- sapply(1:length(to.update.mu), function(x) slideWindow(to.update.mu[x], w_mu[x]) )

    ## Get log prior ratio. Note that the constant parameters will have a prior ratio of 1.
    prop.root.prior <- prior[[1]](prop.root)
    pp <- prop.root.prior - cache.chain$curr.root.prior

    ## Get log likelihood ratio.
    prop.root.lik <- logLikPrunningMCMC(cache.data$X, cache.data$k, cache.data$p, cache.data$nodes, cache.data$des, cache.data$anc, cache.data$mapped.edge
                              , R=cache.chain$chain[[4]], mu=as.vector(prop.root) )
    ll <-  prop.root.lik - cache.chain$lik
    ## Get ratio in log space.
    r <- ll + pp

    ## Acceptance step.
    ## This here need a trick on the for loop. The vcv block is the same as the nex gen.
    if(exp(r) > stats::runif(1)){ ## Accept.
        cat( paste("1; 0; 0; 1; 1; ", prop.root.lik, "\n", sep="") , sep="", file=files[[2]], append=TRUE) ## Always the phylo 1 in this case.
        cache.chain$chain[[1]] <- prop.root
        cache.chain$curr.root.prior <- prop.root.prior
        cache.chain$lik <- prop.root.lik
    } else{                ## Reject.
        cat( paste("0; 0; 0; 1; 1; ", cache.chain$lik, "\n", sep="") , sep="", file=files[[2]], append=TRUE) ## Always the phylo 1 in this case.
    }

    ## Return cache:
    return(cache.chain)
}
