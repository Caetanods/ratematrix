##' Deviance information criterion DIC for mcmcsamples.
##'
##' Function works with the samples from the R matrix MCMC estimate functions.
##' Future implementations will have the DIC calculated alongside the MCMC chain.\cr
##' \cr
##' Issue: DIC is prone to missbehave and reject the model with two or more matrixes even when there is
##'      clearly a better estimate of the parameters. We recommend to use DIC with caution. One can instead
##'      interpret the posterior distribution of matrices and check whether the matrices significantly
##'      overlap.
##' @title Deviance Information Criteria
##' @param out The object with the MCMC info. The return from the MCMC functions.
##' @param post The posterior chain read by the 'read.multi.R.iwish' or 'read.single.R.iwish'.
##' @return numeric. Return the DIC value.
##' @seealso \code{\link{multi.R.iwish.mcmc.R,single.R.iwish.mcmc.R}}
##' @export
dic.mcmc <- function(out, post){
    ## Compute deviance information criterion DIC from mcmcsamples.
    ## Function works with the samples from the R matrix MCMC estimate functions.
    ## The MCMC functions already store the log lik for the parameters. Just need to get them
    ##         and make the calculations. However, still need the original data to make the
    ##         estimate of the likelihood of the posterior mean.
    ## Note that the new version of the readr function will already take out the burn in and
    ##         make a thinning. There is no need to make it here. Thus, 'out$gen' will not have
    ##         the same number of generation of 'post'. Need to use 'post' as it is.
    ## Gelman, A., J. Hwang, and A. Vehtari. 2013. Understanding predictive information
    ##         criteria for Bayesian models. Stat Comput. 24:997–1016.
    ## Future implementation will have DIC calculated alongside the MCMC chain estimate.
    ## ISSUE: I do not recommend you use this function to test hypothesis of one versus two or more
    ##       rate matrices. Simulations showed that DIC is not actuallying showing results congruents to
    ##       reasonable hypothesis testing. A better way to test such and ther hypothesis is on the way.
    ## out = The object with the MCMC info. The return from the MCMC functions.
    ## post = The posterior chain read by the 'read.multi.R.iwish' or 'read.single.R.iwish'.

    ## Get additional objects to calculate the likelihood.
    ## See first lines of 'single.R.iwish.mcmc' for what this is doing.    
    C <- vcv(out$phy)
    X <- out$data[rownames(C),]
    traits <- colnames(out$data)
    y <- matrix( c( as.matrix(out$data) ) )
    n <- length(out$phy$tip.label)
    r <- out$k
    root.post <- post$root
    ss <- dim(root.post)[1]
    b.post <- lapply(1:length(ss), function(y) matrix( sapply(as.vector(root.post[y,]), function(x) rep(x, n) ) ) )

    ## Check the number of matrices fitted in the model.
    if(out$p > 1){

        ## Get the list of C matrices for the multiR.loglik function:
        C.m <- multiC(out$phy) ## phy vcv matrix for multiple R matrices.
        names(C.m) <- colnames(out$phy$mapped.edge) ## Give names to each matrix.
        
        ## Use the multiR.loglik
        ## Number of matrices in the model:
        p <- out$p
        
        ## Get the distribution of parameter estimates from the posterior:
        R.post <- post$matrix

        ## Calculates the mean parameter estimates from the posterior.
        ## This will be used to calculate the first term of equation (7).
        R.mean <- list()
        for(i in 1:p){
            R.mean[[i]] <- Reduce('+', R.post[[i]])/length(R.post[[i]])
        }
        root.mean <- colMeans(root.post)
        b.mean <- matrix( sapply(as.vector(root.mean), function(x) rep(x, n) ) )   
        lik.mean.par <- multiR.loglik(R.m=R.mean, C.m=C.m, y=y, b=b.mean, n=n, r=r, p=p)
        
        ## ## Calculates the mean likelihood across the posterior of parameter estimates.
        ## ## This is already in the output of the MCMC functions.
        ## ## This will be the second term of equation (8).
        lik.mean.post.par <- mean(post$log.lik)
        ## lik.mean.post.par <- mean(post$log.lik[ss])

        ## ## Calculates the p_DIC, effective number of parameters, equation (8):
        ## p.dic <- 2 * ( lik.mean.par - lik.mean.post.par )

        ## ## Calculates the DIC, equation (10.5)
        ## DIC <- (-2 * lik.mean.par) + (2 * p.dic)

        ## Different way, with lik instead of log(lik). The calculation should be the same as before.
        D.bar <- -2 * lik.mean.post.par
        D.hat <- -2 * lik.mean.par
        pD <- D.bar - D.hat
        DIC <- pD + D.bar
        return( unname(DIC) )
        
    } else{        
        ## Get the distribution of parameter estimates from the posterior:
        R.post <- post$matrix
        
        ## Calculates the mean parameter estimates from the posterior.
        ## This will be used to calculate the first term of equation (7).
        R.mean <- Reduce('+', R.post)/length(R.post)
        root.mean <- colMeans(root.post)
        b.mean <- matrix( sapply(as.vector(root.mean), function(x) rep(x, n) ) )
        lik.mean.par <- singleR.loglik(R=R.mean, C=C, y=y, b=b.mean, n=n, r=r)
        
        ## Calculates the mean likelihood across the posterior of parameter estimates.
        ## This will be the second term of equation (8).
        lik.mean.post.par <- mean(post$log.lik)
        
        ## ## Calculates the p_DIC, effective number of parameters, equation (8):
        ## p.dic <- 2 * ( lik.mean.par - lik.mean.post.par )

        ## ## Calculates the DIC, equation (10.5)
        ## DIC <- (-2 * lik.mean.par) + (2 * p.dic)

        ## Different way, with lik instead of log(lik). The calculation should be the same as before.
        D.bar <- -2 * lik.mean.post.par
        D.hat <- -2 * lik.mean.par
        pD <- D.bar - D.hat
        DIC <- pD + D.bar
        return( unname(DIC) )
    }
}
