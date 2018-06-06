##' Take samples from the prior distribution.
##'
##' The prior samples from this function can be used to start the MCMC sampler. See the examples below. \cr
##' \cr
##' If 'sample.sd' is set to FALSE the samples from the standard deviations will be derived from the covariance matrices. If 'sample.sd' is set to TRUE (default) then standard deviations are independently sampled from their own prior distribution and are not derived from the samples of the correlation matrix. Option 'sample.sd = TRUE' is the one used during the MCMC. \cr
##'  \cr
##' The option 'rebuild.R' controls if the samples from the posterior distribution should return the standard deviation separated from the correlation matrix or if these elements should be used to rebuild the covariance matrix. Set 'rebuild.R' to TRUE if you want to obtain the covariance matrices. Otherwise, the 'plotPrior' function works better when 'rebuild.R' is set to FALSE.
##' @title Take samples from the prior distribution
##' @param n number of samples to be generated.
##' @param prior the object with the prior function. See 'makePrior' for more information.
##' @param sample.sd whether the function should sample the vector of standard deviations independently from the correlation matrices. See 'Details'.
##' @param rebuild.R whether the prior sample should return an evolutionary rate matrix rather than a correlation matrix and a vector of standard deviations (default is FALSE). See 'Details'.
##' @return A list with samples from the prior distribution. The structure of this list is the same as required by the parameter 'start' of the 'ratematrixMCMC'.
##' @importFrom corpcor decompose.cov
##' @author Daniel S. Caetano and Luke J. Harmon
##' @export
##' @examples
##' \donttest{
##' data( centrarchidae )
##' dt.range <- t( apply( centrarchidae$data, 2, range ) )
##' ## The step size for the root value can be set given the range we need to sample from:
##' w_mu <- ( dt.range[,2] - dt.range[,1] ) / 10
##' par.sd <- cbind(c(0,0), sqrt( c(10,10) ))
##' prior <- makePrior(r=2, p=2, den.mu="unif", par.mu=dt.range, den.sd="unif", par.sd=par.sd)
##' prior.samples <- samplePrior(n = 1000, prior = prior)
##' start.point <- samplePrior(n=1, prior=prior)
##' ## Plot the prior. Red line shows the sample from the prior that will set the starting 
##' ##        point for the MCMC.
##' plotRatematrix(prior.samples, point.matrix = start.point$matrix, point.color = "red"
##'                , point.wd = 2)
##' handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, prior=prior
##'                          , gen=10000, w_mu=w_mu, dir=tempdir())
##' posterior <- readMCMC(handle, burn = 0.2, thin = 10)
##' ## Again, here the red line shows the starting point of the MCMC.
##' plotRatematrix( posterior, point.matrix = start.point$matrix, point.color = "red"
##'                , point.wd = 2)
##' }
samplePrior <- function(n, prior, sample.sd=TRUE, rebuild.R=FALSE){
    pars <- prior$pars

    ## Sample phylogenetic means:
    mu <- matrix(nrow=n, ncol=pars$r)
    if(pars$den.mu == "unif"){
        for(i in 1:pars$r){
            mu[,i] <- stats::runif(n=n, min=pars$par.mu[i,1], max=pars$par.mu[i,2])
        }
    } else{
        for(i in 1:pars$r){
            mu[,i] <- stats::rnorm(n=n, mean=pars$par.mu[i,1], sd=pars$par.mu[i,2])
        }
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
            if( !is.matrix(pars$Sigma) && !is.list(pars$Sigma) ) stop("Error. Check if the parameter 'Sigma' in function 'makePrior' is of class 'matrix' or class 'list'. Check if length of 'Sigma', if a list, is equal to 'p'.")
        }

        if(sample.sd == TRUE){    
            if(pars$den.sd == "unif"){
                if(pars$p == 1){
                    sd <- stats::runif(n=pars$r, min=pars$par.sd[1], max=pars$par.sd[2])
                } else{
                    sd <- list()
                    for(i in 1:pars$p){ sd[[i]] <- stats::runif(n=pars$r, min=pars$par.sd[i,1], max=pars$par.sd[i,2]) }
                }
            } else{
                if(pars$p == 1){
                    sd <- stats::rlnorm(n=pars$r, meanlog=pars$par.sd[1], sdlog=pars$par.sd[2])
                } else{
                    sd <- list()
                    for(i in 1:pars$p){ sd[[i]] <- stats::rlnorm(n=pars$r, meanlog=pars$par.sd[i,1], sdlog=pars$par.sd[i,2]) }
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
            if( !is.matrix(pars$Sigma) && !is.list(pars$Sigma) ) stop("Error. Check if the parameter 'Sigma' in function 'makePrior' is of class 'matrix' or class 'list'. Check if length of 'Sigma', if a list, is equal to 'p'.")
        }

        if(sample.sd == TRUE){    
            if(pars$den.sd == "unif"){
                if(pars$p == 1){
                    sdvec <- stats::runif(n=pars$r * n, min=pars$par.sd[1], max=pars$par.sd[2])
                    sd <- matrix(sdvec, nrow=n, ncol=pars$r)
                } else{
                    sd <- list()
                    for(i in 1:pars$p){
                        sdvec <- stats::runif(n=pars$r * n, min=pars$par.sd[i,1], max=pars$par.sd[i,2])
                        sd[[i]] <- matrix(sdvec, nrow=n, ncol=pars$r)
                    }
                }
            } else{
                if(pars$p == 1){
                    sdvec <- stats::rlnorm(n=pars$r * n, meanlog=pars$par.sd[1], sdlog=pars$par.sd[2])
                    sd <- matrix(sdvec, nrow=n, ncol=pars$r)
                } else{
                    sd <- list()
                    for(i in 1:pars$p){
                        sdvec <- stats::rlnorm(n=pars$r * n, meanlog=pars$par.sd[i,1], sdlog=pars$par.sd[i,2])
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
    
    ## Check if need to return the R matrix instead of the components.
    if(rebuild.R == TRUE){
        if( n == 1 ){
            if(pars$p == 1){
                R <- rebuild.cov(stats::cov2cor(vcv), v=sd)
            } else{
                R <- lapply(1:pars$p, function(x) rebuild.cov(stats::cov2cor(vcv[[x]]), v=(sd[[x]])^2 ) )
            }
        } else{
            if(pars$p == 1){
                R <- lapply(1:n, function(x) rebuild.cov(stats::cov2cor(vcv[[x]]), v=(sd[x,])^2 ) )
            } else{
                R <- list()
                for(i in 1:pars$p){
                    R[[i]] <- list()
                    for(j in 1:n){
                        R[[i]][[j]] <- rebuild.cov(stats::cov2cor(vcv[[i]][[j]]), v=(sd[[i]][j,])^2 )
                    }
                }
            }
        }
        
        out <- list( root=mu, R=R )
        class( out ) <- "ratematrix_prior_sample_rebuild"
        return( out )
        
    }
    
    out <- list( root=mu, matrix=vcv, sd=sd )
    class( out ) <- "ratematrix_prior_sample"
    return( out )
}
