##' Generates prior densities for the MCMC sampler. Independent priors are defined for the phylogenetic mean, the vector of standard deviations and the structure of correlation. Thus allowing for a wide range of configurations. Priors for the phylogenetic mean and the standard deviations can be uniform of normal (lognormal in the case of the standard deviations). Prior on the matrix of correlations is distributed as a inverse Wishart and can be set to a marginally uniform prior or to be centered is some custom matrix.
##'
##' For a model with two or more rate matrices fitted to the phylogeny different priors can be set for each R matrix. To do this first create a list of 'p' priors created with this function and pass the list to the MCMC sampler. The MCMC will set the prior for the phylogenetic mean (den.mu) as the prior in the first element of the list. The vector of standard deviations and the prior on the correlation will be defined as the other elements of the list in the same order as the R matrices fitted to the data.
##' @title Prior for the MCMC sampler.
##' @param den.mu Set the density type for the phylogenetic mean. One of "unif" or "norm".
##' @param par.mu Numeric vector with length 2. The parameters for the density of phylogenetic means.
##' @param den.sd Set the density type for the vector of standard deviations. One of "unif" or "lnorm".
##' @param par.sd Numeric vector with length 2. The parameters for the density of standard deviations. When set to "unif", then 'par.sd[1]' (the min) need to be a positive value.
##' @param unif.corr Whether the 
##' @param Sigma 
##' @param nu 
##' @return List of density functions to be used as the prior for MCMC functions.
make.prior.zhang <- function(den.mu="unif", par.mu=c(-100,100), den.sd="unif", par.sd=c(0,100), unif.corr=TRUE, Sigma=NULL, nu=NULL){
    ## This should be a much more general prior function than the current implemmented. The point is that for the phylogenetic root and for the vector of standard deviations there are many possibilities for prior distributions. At the moment, two distributions seems good enough. The uniform for both parameters, then the log normal for the vector of standard deviations and the normal for the vector of phylogenetic means.
    ## The parameters for the uniform and the normal are different.
    
    if(!den.mu %in% c("unif","norm") ) stop(" 'den.mu' need to be one of 'unif' or 'norm'.")
    if(!den.sd %in% c("unif","lnorm") ) stop(" 'den.sd' need to be one of 'unif' or 'lnorm'.")
    
    if(den.mu == "unif"){
        mn <- function(x) sum( sapply(x, function(y) dunif(y, min=par.mu[1], max=par.mu[2], log=TRUE)) )
    }
    if(den.mu == "norm"){
        mn <- function(x) sum( sapply(x, function(y) dnorm(y, mean=par.mu[1], sd=par.mu[2], log=TRUE)) )
    }
    
    if(den.sd == "unif"){
        if(par.sd[1] < 0) stop(" 'min' for uniform prior density of standard deviations need to be positive.")
        sd <- function(x) sum( sapply(x, function(y) dunif(y, min=par.sd[1], max=par.sd[2], log=TRUE)) )
    }
    if(den.sd == "lnorm"){
        sd <- function(x) sum( sapply(x, function(y) dlnorm(y, meanlog=par.sd[1], sdlog=par.sd[2], log=TRUE)) )
    }
    
    if(unif.corr == TRUE){
        corr <- function(x) log.diwish(x, v=ncol(x)+1, diag(nrow=ncol(x)) )
    } else{
        if(!is.matrix(Sigma)) stop(" 'Sigma' need to be a matrix. ")
        if(!is.numeric(nu) || nu < ncol(Sigma)) stop(" 'nu' need to be a numeric value =< than the dimension of 'Sigma'. ")
        center <- (nu - ncol(Sigma) -1) * Sigma
        corr <- function(x) log.diwish(x, v=nu, center)
    }
    
    return( list(mn, corr, sd) )
}
