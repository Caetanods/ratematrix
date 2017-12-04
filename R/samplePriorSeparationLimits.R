##' Make samples from the prior density for a given interval.
##'
##' The function will inherit the same parameters from "make.prior.zhang" function.
##'
##' This function is useful to generate a starting value for a MCMC chain and to plot the prior distribution for the model. When generating the starting value make sure to set n=1. To plot the prior distribution you should use a larger sample size to result in a nice-looking plot. For the starting value of a MCMC one can choose whether to sample the starting vector of standard deviations from its prior or to derive from the covariance matrix. Empirical trials show that the MCMC will start easier (the likelihood will behave better) if the vector of standard deviations is computed from the sampled variance covariance matrices rather than sampled independently from the respective prior densities (set sample.sd=FALSE). If using 'sample.sd=TRUE' the MCMC might have dificulties to accept moves related to the evolutionary rate matrices. However, both options are correct for starting the MCMC, it is just a practical matter.
##' @title Sample from the prior
##' @param n Numeric. Number of samples from the prior. Will not work if n == 1.
##' @param prior List. The output from the "make.prior.zhang" function.
##' @param mean.limits a vector with the lower and upper limits for the samples from the prior for the root value. No default value.
##' @param rate.limits a vector with the lower and upper limits for the samples from the prior for the rate matrix. No default value.
##' @param sample.sd Whether the function should sample the vector of standard deviations independently from the samples of the covariance matrices. If set to FALSE then the standard deviations will be derived from the covariance matrices. If set to TRUE (default) then standard deviations are independent from the covariance matrices and are sampled from their correspondent prior density. See 'Details'.
##' @return A list with the parameters draw from the prior. The structure is the same to be used as the "start" parameter of the MCMC chain.
##' @importFrom corpcor decompose.cov
##' @importFrom corpcor rebuild.cov
##' @noRd
samplePriorSeparationLimits <- function(n, prior, mean.limits, rate.limits, sample.sd=TRUE){
    ## At the moment the function works. The problem is that the likelihood to sample a value is low on the regions of the posterior of the MCMC. Thus the rejection sampling takes a lot of time to run. Will set this aside for another time.

    pars <- prior$pars

    ## Sample phylogenetic means:
    mu <- matrix(nrow=n, ncol=pars$r)
    if(pars$den.mu == "unif"){
        for(i in 1:pars$r){ mu[,i] <- stats::runif(n=n, min=mean.limits[1], max=mean.limits[2]) }
    } else{
        for(i in 1:pars$r){
            for(j in 1:n){
                repeat{
                    mu[j,i] <- stats::rnorm(n=1, mean=pars$par.mu[i,1], sd=pars$par.mu[i,2])
                    if( mu[j,i] > mean.limits[1] && mu[j,i] < mean.limits[2] ) break
                }
            }
        }
    }

    ## Sample evolutionary rate matrices from the prior.
    ## The sampling scheme need to be different here. I need to sample the correlation matrix, sample the standard deviation, rebuild the vcv and then check whether the ratematrix is within the limits. So the order of the opperations will change and this version will of course take much more time to run.

    ## Need to define two functions for this. One for the vcv samples and other for the sd samples.
    ## Each of the functions need to return a single vcv matrix or a single vector of standard variations.
    ## So that the functions can be called multiple times to make the rejection sampling.
    sample_vcv <- function(){
        if(pars$unif.corr == TRUE){
            if( pars$p == 1){
                vcv <- riwish(v=pars$r + 1, S=diag(nrow=pars$r))
            } else{
                vcv <- list()
                for(i in 1:pars$p){ vcv[[i]] <- riwish(v=pars$r + 1, S=diag(nrow=pars$r) ) }
            }
        } else{
            if( is.matrix(pars$Sigma) ){
                vcv <- riwish(v=pars$nu, S=pars$Sigma)
            }
            if( is.list(pars$Sigma) ){
                vcv <- list()
                for(i in 1:pars$p){ vcv[[i]] <- riwish(v=pars$nu[i], S=pars$Sigma[[i]]) }
            }
            if( !is.matrix(pars$Sigma) && !is.list(pars$Sigma) ) stop("Error. Check if the parameter 'Sigma' in function 'make.prior.zhang' is of class 'matrix' or class 'list'. Check if length of 'Sigma', if a list, is equal to 'p'.")
        }
        return(vcv) ## A single matrix or a list of matrices.
    }

    sample_sd <- function(){ ## In this case, if the standard deviation is not sampled, just do not call this function.
        if(pars$den.sd == "unif"){
            if(pars$p == 1){
                sd <- stats::runif(n=pars$r, min=pars$par.sd[1], max=pars$par.sd[2])
            } else{
                sd <- list()
                for(i in 1:pars$p){
                    sd[[i]] <- stats::runif(n=pars$r, min=pars$par.sd[i,1], max=pars$par.sd[i,2])
                }
            }
        } else{
            if(pars$p == 1){
                sd <- stats::rlnorm(n=pars$r, meanlog=pars$par.sd[1], sdlog=pars$par.sd[2])
            } else{
                sd <- list()
                for(i in 1:pars$p){
                    sd[[i]] <- stats::rlnorm(n=pars$r, meanlog=pars$par.sd[i,1], sdlog=pars$par.sd[i,2])
                }
            }
        }
        return(sd) ## A single vector or a list of vectors.
    }

    check_limit <- function(X){
        ## To check if every cell of the matrix is within the limits.
        mm <- c( X[ upper.tri(X, diag=TRUE) ] )
        low <- mm > rate.limits[1]
        high <- mm < rate.limits[2]
        if( sum(low) == length(mm) && sum(high) == length(mm) ){
            return(TRUE)
        } else{
            return(FALSE)
        }
    }

    ## Now need to perform the rejection sampling:
    if( sample.sd == TRUE ){
        prior.samples <- list()
        for( i in 1:n ){
            repeat{
                vcv <- sample_vcv()
                sd <- sample_sd()
                corr <- decompose.cov(vcv)$r
                vv <- (sd)^2
                prior.samples[[i]] <- rebuild.cov(r=corr, v=vv) ## Reject
                if( check_limit( prior.samples[[i]] ) ) break
            }
        }
    } else{
        prior.samples <- list()
        for( i in 1:n ){
            repeat{
                prior.samples[[i]] <- sample_vcv()
                if( check_limit( prior.samples[[i]] ) ) break
            }
        }
    }

    out <- list( mu=mu, matrix=prior.samples )
    return( out ) ## This output is different from the other type of prior samples. Might need to adapt.
}
