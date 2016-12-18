##' Generates samples from the prior distribution used in the MCMC chain.
##'
##' The function will inherit the same parameters from "make.prior.zhang" function.
##'
##' 
##' @title Take samples from the prior distribution
##' @param n number of samples to be generated.
##' @param prior the object with the prior function. See 'makePrior' for more information.
##' @param sample.sd whether the function should sample the vector of standard deviations independently from the correlation matrices. If set to FALSE, the samples for standard deviation will be derived from the covariance matrices. If set to TRUE (default) then standard deviations are independently sampled from the own prior distribution and are not derived from the samples of the correlation matrix.
##' @return A list with samples from the prior distribution. The structure of this list is the same as required by the parameter 'start' of the 'ratematrixMCMC'.
##' @importFrom corpcor decompose.cov
##' @author Daniel S. Caetano and Luke J. Harmon
##' @export
##' @seealso \code{\link{ plotRatematrix }} and \code{\link{ plotRootValue }} for plotting the samples from the prior.
##' @examples
##' \donttest{
##' par.mu <- rbind( c(-10, 10), c(-10, 10), c(-10, 10) )
##' par.sd <- c(0, 10)
##' prior <- makePrior(r=3, p=1, par.mu=par.mu, par.sd=par.sd)
##' ## Sample prior. Standard deviation sampled independently.
##' samples.ind <- samplePrior(n=1000, prior=prior, sample.sd=TRUE)
##' ## Now plot.
##' plotRatematrix(samples.ind, set.xlim=c(-100,100))
##' ## Sample prior. Now standard deviation is NOT sampled, just derived from the vcv matrices
##' ##    sampled from the inverse-Wishart distribution.
##' samples.cor <- samplePrior(n=1000, prior=prior, sample.sd=FALSE)
##' ## Note how the inverse-Wishart is much more informative than the 'separation strategy'.
##' ##    The x axis is the same.
##' plotRatematrix(samples.cor, set.xlim=c(-100,100))
##' 
##' ## We can also use the sample from the prior to start the MCMC:
##' par.mu <- rbind( c(-10, 10), c(-10, 10) )
##' par.sd <- rbind( c(0, 10), c(0, 10) )
##' prior <- makePrior(r=2, p=2, par.mu=par.mu, par.sd=par.sd)
##' ## Need to take a single sample to be used as the starting point.
##' start <- samplePrior(n=1, prior=prior)
##' data(centrarchidae)
##' handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=1000, start=start)
##' }
samplePrior <- function(n, prior, sample.sd=TRUE){
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
    out <- list( mu=mu, matrix=vcv, sd=sd )
    class( out ) <- "ratematrix_prior_sample"
    return( out )    
}
