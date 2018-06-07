##' Description!
##'
##' Details!
##' @title Sigma step using Zhang strategy
##' @param cache.data 
##' @param cache.chain 
##' @param prior 
##' @param v 
##' @param w_sd 
##' @param w_mu 
##' @param iter 
##' @param count 
##' @param n.phy 
##' @return The chain cache.
##' @author daniel
##' @importFrom corpcor decompose.cov rebuild.cov
##' @noRd
makePropMultSigmaList <- function(cache.data, cache.chain, prior, v, w_sd, w_mu, files, n.phy) {
    ## This is going to be the step for the correlation matrix and the vector of standard deviations.
    ## The moves for the correlation matrix now are independent of the moves for the standard deviations.
    ## Thus, I need to sample which move to make at each call of the function.
    ## In a second implementation, would be better to separate this function into two independent update functions maybe.
    
    which.phy <- sample(1:n.phy, size = 1) ## The index for the phy used for the log.lik.
    up <- sample(1:2, size=1) ## Updates the standard deviations or the correlation structure.
    ## Random draw one of the p R matrices to update:
    Rp <- sample(1:cache.data$p, size=1)

    if( up == 1 ){
        
        ## Update the vector of standard deviations.
        to.update.sd <- cache.chain$chain[[3]][[Rp]]
        prop.sd <- sapply(1:length(to.update.sd), function(x) slideWindowPositive(to.update.sd[x], w_sd[x]) )
        ## Need to put the parameters together for the prior.
        prop.sd.full <- cache.chain$chain[[3]]
        prop.sd.full[[Rp]] <- prop.sd
        prop.sd.prior <- prior[[3]](prop.sd.full) ## The third prior function.
        pp <- prop.sd.prior - cache.chain$curr.sd.prior ## Two numeric quantities.

        ## Rebuild the matrix to calculate the likelihood.
        ## No need for the Jacobian in this move.
        decom <- decompose.cov( cache.chain$chain[[4]][[Rp]] )
        ## prop.vcv is the covariance matrix to calculate the likelihood and the parameter for the posterior.
        prop.vcv <- rebuild.cov( r=decom$r, v=prop.sd^2 ) ## We are making moves to the standard deviation.
        ## This part will not work with the 'log.dmvnorm' function from the 'ratematrix' package.
        ## Only work with the one defined here.
        Rlik <- cache.chain$chain[[4]] ## This log lik need both matrices.
        Rlik[[Rp]] <- prop.vcv
        prop.sd.lik <- logLikPrunningMCMC(cache.data$X, cache.data$k, cache.data$p, cache.data$nodes[[which.phy]], cache.data$des[[which.phy]]
                                        , cache.data$anc[[which.phy]], cache.data$mapped.edge[[which.phy]], R=Rlik
                                        , mu=as.vector( cache.chain$chain[[1]] ) )
        ## prop.sd.lik <- logDensityMvNorm(cache.data$X, mu=cache.chain$chain[[iter-1]][[1]], sigma=prop.vcv)
        ll <-  prop.sd.lik - cache.chain$lik
        
        ## Get ratio in log space.
        r <- ll + pp

        ## Acceptance step.
        ## This here need a trick on the for loop. The vcv block is the same as the nex gen.
        if(exp(r) > stats::runif(1)){ ## Accept.
            cat( paste("1; 0; ", Rp, "; 0; ", which.phy, "; ", prop.sd.lik, "\n", sep="") , sep="", file=files[[2]], append=TRUE) ## Rp = the regime updated.
            cache.chain$chain[[3]][[Rp]] <- prop.sd
            cache.chain$chain[[4]][[Rp]] <- prop.vcv
            cache.chain$curr.sd.prior <- prop.sd.prior
            cache.chain$lik <- prop.sd.lik
        } else{                ## Reject.
            cat( paste("0; 0; ", Rp, "; 0; ", which.phy, "; ", cache.chain$lik, "\n", sep="") , sep="", file=files[[2]], append=TRUE) ## Rp = the regime updated.
        }
        
    }

    if( up == 2 ){
        ## Update the correlation matrix.
        prop.r <- makePropIWish(cache.chain$chain[[2]][[Rp]], k=cache.data$k, v=v)
        prop.r.full <- cache.chain$chain[[2]]
        prop.r.full[[Rp]] <- prop.r
        prop.r.prior <- prior[[2]](prop.r.full) ## The second prior function. On the expanded parameters, the covariance matrix.
        pp <- prop.r.prior - cache.chain$curr.r.prior

        ## Rebuild the matrix to calculate the likelihood:
        decom <- decompose.cov( prop.r )
        ## We use the independent vector of standard deviations to calculate the likelihood.
        ## prop.vcv is the covariance matrix to calculate the likelihood and the parameter for the posterior.
        prop.vcv <- rebuild.cov( r=decom$r, v=cache.chain$chain[[3]][[Rp]]^2 )
        Rlik <- cache.chain$chain[[4]] ## This log lik need both matrices.
        Rlik[[Rp]] <- prop.vcv
        prop.r.lik <- logLikPrunningMCMC(cache.data$X, cache.data$k, cache.data$p, cache.data$nodes[[which.phy]], cache.data$des[[which.phy]]
                                       , cache.data$anc[[which.phy]], cache.data$mapped.edge[[which.phy]], R=Rlik
                                       , mu=as.vector( cache.chain$chain[[1]] ) )
        ll <- prop.r.lik - cache.chain$lik
        ## The hastings ratio.
        hh <- hastingsDensity(curr.vcv=cache.chain$chain[[2]][[Rp]], prop.vcv=prop.r, p=cache.data$k, v=v)
        ## Need the Jacobian, this is in function of the vector of VARIANCES, not the standard deviations.
        prop.r.jacobian <- sum( sapply(1:cache.data$k, function(x) log( decom$v[x]) ) ) * log( (cache.data$k-1)/2 )
        jj <- prop.r.jacobian - cache.chain$curr.r.jacobian[[Rp]]

        ## Get ratio in log space. Log lik, log prior, log hastings and log jacobian.
        r <- ll + pp + hh + jj

        ## Acceptance step.
        ## This here need a trick on the for loop. The vcv block is the same as the nex gen.
        if(exp(r) > stats::runif(1)){ ## Accept.
            cat( paste("1; ", Rp, "; 0; 0; ", which.phy, "; ", prop.r.lik, "\n", sep="") , sep="", file=files[[2]], append=TRUE) ## Rp = the regime updated.
            cache.chain$chain[[2]][[Rp]] <- prop.r
            cache.chain$chain[[4]][[Rp]] <- prop.vcv
            cache.chain$curr.r.prior <- prop.r.prior
            cache.chain$curr.r.jacobian[[Rp]] <- prop.r.jacobian
            cache.chain$lik <- prop.r.lik
        } else{                ## Reject.
            cat( paste("0; ", Rp, "; 0; 0; ", which.phy, "; ", cache.chain$lik, "\n", sep="") , sep="", file=files[[2]], append=TRUE) ## Rp = the regime updated.
        }
        
    }

## Return cache:
return(cache.chain)
}
