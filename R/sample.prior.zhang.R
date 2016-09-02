##' Make samples from the prior density.
##'
##' The function will inherit the same parameters from "make.prior.zhang" function.
##'
##' This function is useful to generate a starting value for a MCMC chain and to plot the prior distribution for the model. When generating the starting value make sure to set n=1. To plot the prior distribution you should use a larger sample size to result in a nice-looking plot. For the starting value of a MCMC one can choose whether to sample the starting vector of standard deviations from its prior or to derive from the covariance matrix. Empirical trials show that the MCMC will start easier (the likelihood will behave better) if the vector of standard deviations is computed from the sampled variance covariance matrices rather than sampled independently from the respective prior densities (set sample.sd=FALSE). If using 'sample.sd=TRUE' the MCMC might have dificulties to accept moves related to the evolutionary rate matrices. However, both options are correct for starting the MCMC, it is just a practical matter.
##' @title Sample from the prior
##' @param n Numeric. Number of samples from the prior.
##' @param prior List. The output from the "make.prior.zhang" function.
##' @param sample.sd Whether the function should sample the vector of standard deviations independently from the samples of the covariance matrices. If set to FALSE then the standard deviations will be derived from the covariance matrices. If set to TRUE (default) then standard deviations are independent from the covariance matrices and are sampled from their correspondent prior density. See 'Details'.
##' @return A list with the parameters draw from the prior. The structure is the same to be used as the "start" parameter of the MCMC chain.
##' @importFrom corpcor decompose.cov
##' @export
sample.prior.zhang <- function(n, prior, sample.sd=TRUE){
    pars <- prior$pars

    ## Sample phylogenetic means:
    mu <- matrix(nrow=n, ncol=pars$r)
    if(pars$den.mu == "unif"){
        for(i in 1:pars$r){ mu[,i] <- runif(n=n, min=pars$par.mu[i,1], max=pars$par.mu[i,2]) }
    } else{
        for(i in 1:pars$r){ mu[,i] <- rnorm(n=n, mean=pars$par.mu[i,1], sd=pars$par.mu[i,2]) }
    }

    if(n == 1){ mu <- as.numeric(mu[1,]) }

    ## Sample evolutionary rate matrices. They are composed by a vector of standard deviations and a covariance matrix.
    ## Standard deviations come from a separate and independent prior. So they need to be sampled indepedently too.
    ## Would be cool to add an option to only sample the evolutionary rate matrix and derive sd from those matrices.

    if(n == 1){ ## A single sample. Good for the start of the MCMC.
        if(pars$unif.corr == TRUE){
            if( pars$p == 1){
                vcv <- riwish(v=pars$r + 1, S=diag(nrow=pars$r) )
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

        if(sample.sd == TRUE){    
            if(pars$den.sd == "unif"){
                if(pars$p == 1){
                    sd <- runif(n=pars$r, min=pars$par.sd[1], max=pars$par.sd[2])
                } else{
                    sd <- list()
                    for(i in 1:pars$p){ sd[[i]] <- runif(n=pars$r, min=pars$par.sd[i,1], max=pars$par.sd[i,2]) }
                }
            } else{
                if(pars$p == 1){
                    sd <- rlnorm(n=pars$r, meanlog=pars$par.sd[1], sdlog=pars$par.sd[2])
                } else{
                    sd <- list()
                    for(i in 1:pars$p){ sd[[i]] <- rlnorm(n=pars$r, meanlog=pars$par.sd[i,1], sdlog=pars$par.sd[i,2]) }
                }
            }
        } else{
            if(pars$p == 1){
                sd <- sqrt( decompose.cov(vcv)$v )
            } else{
                sd <- list()
                for(i in 1:pars$p){ sd[[i]] <- sqrt( decompose.cov(vcv[[i]])$v ) }
            }
        }
    }
    if(n > 1){ ## Several samples. Good for plotting.
        if(pars$unif.corr == TRUE){
            if( pars$p == 1){
                vcv <- lapply(1:n, function(x) riwish(v=pars$r + 1, S=diag(nrow=pars$r) ) )
            } else{
                vcv <- list()
                for(i in 1:pars$p){ vcv[[i]] <- lapply(1:n, function(x) riwish(v=pars$r + 1, S=diag(nrow=pars$r) ) ) }
            }
        } else{
            if( is.matrix(pars$Sigma) ){
                vcv <- lapply(1:n, function(x) riwish(v=pars$nu, S=pars$Sigma) )
            }
            if( is.list(pars$Sigma) ){
                vcv <- list()
                for(i in 1:pars$p){ vcv[[i]] <- lapply(1:n, function(x) riwish(v=pars$nu[i], S=pars$Sigma[[i]]) ) }
            }
            if( !is.matrix(pars$Sigma) && !is.list(pars$Sigma) ) stop("Error. Check if the parameter 'Sigma' in function 'make.prior.zhang' is of class 'matrix' or class 'list'. Check if length of 'Sigma', if a list, is equal to 'p'.")
        }

        if(sample.sd == TRUE){    
            if(pars$den.sd == "unif"){
                if(pars$p == 1){
                    sdvec <- runif(n=pars$r * n, min=pars$par.sd[1], max=pars$par.sd[2])
                    sd <- matrix(sdvec, nrow=n, ncol=pars$r)
                } else{
                    sd <- list()
                    for(i in 1:pars$p){
                        sdvec <- runif(n=pars$r * n, min=pars$par.sd[i,1], max=pars$par.sd[i,2])
                        sd[[i]] <- matrix(sdvec, nrow=n, ncol=pars$r)
                    }
                }
            } else{
                if(pars$p == 1){
                    sdvec <- rlnorm(n=pars$r * n, meanlog=pars$par.sd[1], sdlog=pars$par.sd[2])
                    sd <- matrix(sdvec, nrow=n, ncol=pars$r)
                } else{
                    sd <- list()
                    for(i in 1:pars$p){
                        sdvec <- rlnorm(n=pars$r * n, meanlog=pars$par.sd[i,1], sdlog=pars$par.sd[i,2])
                        sd[[i]] <- matrix(sdvec, nrow=n, ncol=pars$r)
                    }
                }
            }
        } else{
            if(pars$p == 1){
                sd <- t( sapply(1:n, function(x) sqrt( decompose.cov(vcv[[x]])$v ) ) )
            } else{
                sd <- list()
                for(i in 1:pars$p){
                    sd[[i]] <- t( sapply(1:n, function(x) sqrt( decompose.cov(vcv[[i]][[x]])$v ) ) )
                }
            }
        }
    }
    return( list( mu=mu, vcv=vcv, sd=sd ) )    
}
