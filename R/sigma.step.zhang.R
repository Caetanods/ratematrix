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
sigma.step.zhang <- function(cache.data, cache.chain, prior, w_sd, w_mu, v, iter, count) {
    ## This is going to be the step for the correlation matrix and the vector of standard deviations.
    ## The moves for the correlation matrix now are independent of the moves for the standard deviations.
    ## Thus, I need to sample which move to make at each call of the function.
    ## In a second implementation, would be better to separate this function into two independent update functions maybe.

    up <- sample(1:2, size=1)

    if( up == 1 ){
        ## Update the vector of standard deviations.
        prop.sd <- sapply(cache.chain$chain[[iter-1]][[3]], function(x) sliding.window.positive(x, w_sd) )
        prop.sd.prior <- prior[[3]]( prop.sd ) ## The third prior function. New prior works on list format.
        pp <- prop.sd.prior - cache.chain$curr.sd.prior

        ## Rebuild the matrix to calculate the likelihood.
        ## No need for the Jacobian in this move.
        decom <- decompose.cov( cache.chain$chain[[iter-1]][[2]] )
        ## prop.vcv is the covariance matrix to calculate the likelihood and the parameter for the posterior.
        prop.vcv <- rebuild.cov( r=decom$r, v=prop.sd^2 ) ## We are making moves to the standard deviation.
        ## This part will not work with the 'log.dmvnorm' function from the 'ratematrix' package.
        ## Only work with the one defined here.
        prop.sd.lik <- singleR.loglik(data=cache.data, chain=cache.chain, root=as.vector( cache.chain$chain[[iter-1]][[1]] )
                                    , R=prop.vcv)
        ## prop.sd.lik <- log.dmvnorm(cache.data$X, mu=cache.chain$chain[[iter-1]][[1]], sigma=prop.vcv)
        ll <-  prop.sd.lik - cache.chain$lik[iter-1]
        
        ## Get ratio in log space.
        r <- ll + pp

        ## Acceptance step.
        ## This here need a trick on the for loop. The vcv block is the same as the nex gen.
        if(exp(r) > runif(1)){ ## Accept.
            cache.chain$chain[[iter]] <- cache.chain$chain[[iter-1]]
            cache.chain$chain[[iter]][[3]] <- prop.sd
            cache.chain$chain[[iter]][[4]] <- prop.vcv
            cache.chain$curr.sd.prior <- prop.sd.prior
            cache.chain$acc[count] <- 3
            cache.chain$lik[iter] <- prop.sd.lik
        } else{                ## Reject.
            cache.chain$chain[[iter]] <- cache.chain$chain[[iter-1]]
            cache.chain$acc[count] <- 0
            cache.chain$lik[iter] <- cache.chain$lik[iter-1]
        }
        
    }

    if( up == 2 ){
        ## Update the correlation matrix.
        prop.r <- make.prop.iwish(cache.chain$chain[[iter-1]][[2]], k=cache.data$k, v=v)
        prop.r.prior <- prior[[2]]( prop.r ) ## The second prior function. On the expanded parameters, the covariance matrix.
        pp <- prop.r.prior - cache.chain$curr.r.prior

        ## Rebuild the matrix to calculate the likelihood:
        decom <- decompose.cov( prop.r )
        ## We use the independent vector of standard deviations to calculate the likelihood.
        ## prop.vcv is the covariance matrix to calculate the likelihood and the parameter for the posterior.
        prop.vcv <- rebuild.cov( r=decom$r, v=cache.chain$chain[[iter-1]][[3]]^2 )
        prop.r.lik <- singleR.loglik(data=cache.data, chain=cache.chain, root=as.vector( cache.chain$chain[[iter-1]][[1]] )
                                   , R=prop.vcv)
        ll <- prop.r.lik - cache.chain$lik[iter-1]
        ## The hastings ratio.
        hh <- h.density(curr.vcv=cache.chain$chain[[iter-1]][[2]], prop.vcv=prop.r, p=cache.data$k, v=v)
        ## Need the Jacobian, this is in function of the vector of VARIANCES, not the standard deviations.
        prop.r.jacobian <- sum( sapply(1:cache.data$k, function(x) log( decom$v[x]) ) ) * log( (cache.data$k-1)/2 )
        jj <- prop.r.jacobian - cache.chain$curr.r.jacobian

        ## Get ratio in log space. Log lik, log prior, log hastings and log jacobian.
        r <- ll + pp + hh + jj

        ## Acceptance step.
        ## This here need a trick on the for loop. The vcv block is the same as the nex gen.
        if(exp(r) > runif(1)){ ## Accept.
            cache.chain$chain[[iter]] <- cache.chain$chain[[iter-1]]
            cache.chain$chain[[iter]][[2]] <- prop.r
            cache.chain$chain[[iter]][[4]] <- prop.vcv
            cache.chain$curr.r.prior <- prop.r.prior
            cache.chain$curr.r.jacobian <- prop.r.jacobian
            cache.chain$acc[count] <- 2
            cache.chain$lik[iter] <- prop.r.lik
        } else{                ## Reject.
            cache.chain$chain[[iter]] <- cache.chain$chain[[iter-1]]
            cache.chain$acc[count] <- 0
            cache.chain$lik[iter] <- cache.chain$lik[iter-1]
        }
        
    }

## Return cache:
return(cache.chain)
}
