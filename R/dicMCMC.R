##' Deviance information criterion (DIC) for MCMC models.
##'
##' Function works with the samples from the R matrix MCMC estimate functions.
##' \cr
##' \cr
##' Important note: DIC uses the mean of the posterior as a summary statistics for the test. If the posterior of parameter estimates is not well summarized by the mean, then the results of the DIC test will have no meaning.
##' @title Deviance Information Criteria
##' @param out The object with the MCMC info. The return from the MCMC functions.
##' @param post The posterior chain read by the 'read.multi.R.iwish' or 'read.single.R.iwish'.
##' @return numeric. Return the DIC value.
##' @seealso \code{\link{multi.R.iwish.mcmc.R single.R.iwish.mcmc.R}}
##' @references Gelman, A., J. Hwang, and A. Vehtari. 2013. Understanding predictive information criteria for Bayesian models. Stat Comput. 24:997–1016.
##' @importFrom ape vcv.phylo as.phylo
##' @importFrom phytools multiC
##' @importFrom corpcor make.positive.definite
##' @export
dicMCMC <- function(out, post){
    cache.data <- list()
    cache.data$C <- vcv.phylo(out$phy)
    cache.data$X <- out$data[rownames(cache.data$C),]
    cache.data$n <- length(out$phy$tip.label)
    cache.data$r <- out$k
    D <- matrix(0, nrow = cache.data$n * cache.data$r, ncol = cache.data$r)
    for(i in 1:cache.data$r) D[((cache.data$n*(i-1))+1):(cache.data$n*i),i] <- 1
    cache.data$D <- D
    traits <- colnames(out$data)
    y <- matrix( c( as.matrix(out$data) ) )
    root.post <- post$root
    ss <- dim(root.post)[1]

    ## Check the number of matrices fitted in the model.
    if(out$p > 1){

        ## Get the list of C matrices for the multiR.loglik function:
        cache.data$C.m <- multiC(out$phy) ## phy vcv matrix for multiple R matrices.
        names(cache.data$C.m) <- colnames(out$phy$mapped.edge) ## Give names to each matrix.
        
        ## Use the multiR.loglik
        ## Number of matrices in the model:
        cache.data$p <- out$p
        
        ## Get the distribution of parameter estimates from the posterior:
        R.post <- post$matrix

        ## Calculates the mean parameter estimates from the posterior.
        ## This will be used to calculate the first term of equation (7).
        R.mean <- list()
        for(i in 1:cache.data$p){
            RR <- Reduce('+', R.post[[i]])/length(R.post[[i]])
            if( det(RR) == 0 ) RR <- make.positive.definite(RR)
            R.mean[[i]] <- RR
        }
        
        root.mean <- colMeans(root.post)
        ## Now the likelihood function depends on the caches. So need to update this.
        lik.mean.par <- logLikMultRegime(data=cache.data, chain=NULL, root=root.mean, R=R.mean)
        
        ## Calculates the mean likelihood across the posterior of parameter estimates.
        ## This is already in the output of the MCMC functions.
        ## This will be the second term of equation (8).
        lik.mean.post.par <- mean(post$log.lik)

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

        ## Need to check if the matrix from the mean parameter estimates is positive-definite.
        ## if not, generate the nearest positive-definite matrix.
        if( det(R.mean) == 0 ) R.mean <- make.positive.definite(R.mean)

        root.mean <- colMeans(root.post)
        cache.data$phy <- as.phylo(out$phy)
        lik.mean.par <- logLikSingleRegime(data=cache.data, chain=NULL, root=root.mean, R=R.mean)
        
        ## Calculates the mean likelihood across the posterior of parameter estimates.
        ## This will be the second term of equation (8).
        lik.mean.post.par <- mean(post$log.lik)
        
        ## Different way, with lik instead of log(lik). The calculation should be the same as before.
        D.bar <- -2 * lik.mean.post.par
        D.hat <- -2 * lik.mean.par
        pD <- D.bar - D.hat
        DIC <- pD + D.bar
        return( unname(DIC) )
    }
}
