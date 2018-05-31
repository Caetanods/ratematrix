##' Description!
##'
##' Details!
##' @title Sigma step using Zhang strategy
##' @param cache.data Description
##' @param cache.chain Description
##' @param prior Description
##' @param w_sd Description
##' @param w_mu Description
##' @param v Description
##' @param iter Description
##' @param count Description
##' @return The chain cache.
##' @importFrom corpcor decompose.cov rebuild.cov
##' @noRd
makePropSingleSigmaList <- function(cache.data, cache.chain, prior, w_sd, w_mu, v, files, phy, n.phy) {

    which.phy <- sample(1:n.phy, size = 1) ## The index for the phy used for the log.lik.

    mat.prob <- cache.data$prop[2:3] / sum(cache.data$prop[2:3])
    up <- sample(1:2, size=1, prob = mat.prob)

    if( up == 1 ){
        ## Update the vector of standard deviations.
        to.update.sd <- cache.chain$chain[[3]]
        ## multi.prop[1,] is the proposed vector.
        ## multi.prop[2,] is the prop ratio vector.
        multi.prop <- sapply(1:length(to.update.sd), function(x) multiplierProposal(x = to.update.sd[x], a = w_sd[x]) )
        prop.sd <- multi.prop[1,]
        prop.sd.prior <- prior[[3]]( prop.sd ) ## The third prior function. New prior works on list format.
        pp <- prop.sd.prior - cache.chain$curr.sd.prior

        ## Rebuild the matrix to calculate the likelihood.
        ## No need for the Jacobian in this move.
        decom <- decompose.cov( cache.chain$chain[[4]] )
        ## prop.vcv is the covariance matrix to calculate the likelihood and the parameter for the posterior.
        prop.vcv <- rebuild.cov( r=decom$r, v=prop.sd^2 ) ## We are making moves to the standard deviation.
        ## This part will not work with the 'log.dmvnorm' function from the 'ratematrix' package.
        ## Only work with the one defined here.
        prop.sd.lik <- logLikSingleRegime(data=cache.data, chain=cache.chain, phy=phy[[which.phy]]
                                        , root=as.vector( cache.chain$chain[[1]] )
                                        , R=prop.vcv)
        ## prop.sd.lik <- logDensityMvNorm(cache.data$X, mu=cache.chain$chain[[iter-1]][[1]], sigma=prop.vcv)
        ll <-  prop.sd.lik - cache.chain$lik
        
        ## Get ratio in log space.
        ## multi.prop is the proposal ratio for the multiplier proposal.
        r <- ll + pp + sum( multi.prop[2,] )

        ## Acceptance step.
        ## This here need a trick on the for loop. The vcv block is the same as the nex gen.
        if(exp(r) > stats::runif(1)){ ## Accept.
            cat( paste("1; 0; 1; 0; ", which.phy, "; ", prop.sd.lik, "\n", sep="") , sep="", file=files[[2]], append=TRUE)
            cache.chain$chain[[3]] <- prop.sd
            cache.chain$chain[[4]] <- prop.vcv
            cache.chain$curr.sd.prior <- prop.sd.prior
            cache.chain$lik <- prop.sd.lik
        } else{                ## Reject.
            cat( paste("0; 0; 1; 0; ", which.phy, "; ", cache.chain$lik, "\n", sep="") , sep="", file=files[[2]], append=TRUE)
        }
        
    }

    if( up == 2 ){
        ## Update the correlation matrix.
        prop.r <- makePropIWish(cache.chain$chain[[2]], k=cache.data$k, v=v)
        prop.r.prior <- prior[[2]]( prop.r ) ## The second prior function. On the expanded parameters, the covariance matrix.
        pp <- prop.r.prior - cache.chain$curr.r.prior

        ## Rebuild the matrix to calculate the likelihood:
        decom <- decompose.cov( prop.r )
        ## We use the independent vector of standard deviations to calculate the likelihood.
        ## prop.vcv is the covariance matrix to calculate the likelihood and the parameter for the posterior.
        prop.vcv <- rebuild.cov( r=decom$r, v=cache.chain$chain[[3]]^2 )
        prop.r.lik <- logLikSingleRegime(data=cache.data, chain=cache.chain, phy=phy[[which.phy]]
                                       , root=as.vector( cache.chain$chain[[1]] )
                                       , R=prop.vcv)
        ll <- prop.r.lik - cache.chain$lik
        ## The hastings ratio.
        hh <- hastingsDensity(curr.vcv=cache.chain$chain[[2]], prop.vcv=prop.r, p=cache.data$k, v=v)
        ## Need the Jacobian, this is in function of the vector of VARIANCES, not the standard deviations.
        prop.r.jacobian <- sum( sapply(1:cache.data$k, function(x) log( decom$v[x]) ) ) * log( (cache.data$k-1)/2 )
        jj <- prop.r.jacobian - cache.chain$curr.r.jacobian

        ## Get ratio in log space. Log lik, log prior, log hastings and log jacobian.
        r <- ll + pp + hh + jj

        ## Acceptance step.
        ## This here need a trick on the for loop. The vcv block is the same as the nex gen.
        if(exp(r) > stats::runif(1)){ ## Accept.
            cat( paste("1; 1; 0; 0; ", which.phy, "; ", prop.r.lik, "\n", sep="") , sep="", file=files[[2]], append=TRUE)
            cache.chain$chain[[2]] <- prop.r
            cache.chain$chain[[4]] <- prop.vcv
            cache.chain$curr.r.prior <- prop.r.prior
            cache.chain$curr.r.jacobian <- prop.r.jacobian
            cache.chain$lik <- prop.r.lik
        } else{                ## Reject.
            cat( paste("0; 1; 0; 0; ", which.phy, "; ", cache.chain$lik, "\n", sep="") , sep="", file=files[[2]], append=TRUE)
        }
        
    }

## Return cache:
return(cache.chain)
}
