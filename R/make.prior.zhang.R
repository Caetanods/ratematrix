##' Generates prior densities for the MCMC sampler. Independent priors are defined for the phylogenetic mean, the vector of standard deviations and the structure of correlation. Thus allowing for a wide range of configurations. Priors for the phylogenetic mean and the standard deviations can be uniform or normal (lognormal in the case of the standard deviations). Prior on the matrix of correlations is distributed as an inverse Wishart and can be set to a marginally uniform prior or to be centered is some custom matrix.
##'
##' For a model with two or more rate matrices fitted to the phylogeny different priors can be set for each R matrix. To do this first create a list of 'p' priors created with this function and pass the list to the MCMC sampler. The MCMC will set the prior for the phylogenetic mean (den.mu) as the prior in the first element of the list. The vector of standard deviations and the prior on the correlation will be defined as the other elements of the list in the same order as the R matrices fitted to the data.
##' @title Prior for the MCMC sampler.
##' @param r Number of traits in the model.
##' @param p Number of evolutionary rate matrices fitted to the phylogeny.
##' @param den.mu Set the density type for the phylogenetic mean. One of "unif" (default) or "norm".
##' @param par.mu The parameters for the prior density on the vector of phylogenetic means. Matrix with number of rows equal 'r'. When the density is set to "unif" then par.mu[,1] is the minimum values and par.mu[,2] is the maximum values for each trait. When the density is set to "norm̉" then par.mu[,1] is the mean values and par.mu[,2] is the standard deviation values for the set of normal densities around the vector of phylogenetic means.
##' @param den.sd Set the density type for the vector of standard deviations. One of "unif" or "lnorm".
##' @param par.sd The parameters for the density of standard deviations. Matrix with number of rows equal to 'p'. When "den.sd" is set to "unif", then 'par.sd[,1]' (the minimum) need to be a vector of positive values and 'par.sd[,2]' is the vector of maximum values. When "den.sd" is set to "lnorm" then 'par.sd[,1]' is the vector of log(means) for the density and 'par.sd[,2]' is the vector of log(sd) for the distributions.
##' @param unif.corr Whether the correlation structure of the prior distribution on the Sigma matrix is flat.
##' @param Sigma The scale matrix for the inverse Wishart distribution used to sample the evolutionary correlation among traits. This is the parameter to center the prior distribution of the correlations in a given value. Need to be a list with length equal to 'p'. When 'p==1' then is needs to be a matrix.
##' @param nu The degrees of freedom parameter of the inverse Wishart distribution. Larger values will make the distribution more tight around the scale matrix. Need to be a vector of length equal to 'p'.
##' @return List of density functions to be used as the prior for MCMC functions.
##' @export
make.prior.zhang <- function(r, p, den.mu="unif", par.mu, den.sd="unif", par.sd, unif.corr=TRUE, Sigma=NULL, nu=NULL){

    ## Check if the distributions are known. This assumes all priors for the same parameter shares a distribution (with different parameters).
    ## Could extend this in the future.
    if(!den.mu %in% c("unif","norm") ) stop(" 'den.mu' need to be one of 'unif' or 'norm'.")
    if(!den.sd %in% c("unif","lnorm") ) stop(" 'den.sd' need to be one of 'unif' or 'lnorm'.")

    ## Check the number of parameters sent to the function. Those need to match.
    lpmu <- nrow(par.mu)
    if(!r == lpmu) stop("number of rows of 'par.mu' need to be equal to the number of traits in the model (r).")
    if( is.vector(par.sd) && p > 1 ) stop("p need to be equal to the regimes fitted to the tree.")
    if( is.matrix(par.sd) && !p == nrow(par.sd) ) stop("p need to be equal to the regimes fitted to the tree. So p == nrow(par.sd).")    
    ## lpsd <- nrow(par.sd)
    ## if(!p == lpsd) stop("number of rows of 'par.sd' need to be equal to the number of R matrices fitted to the data (p).")
    if( unif.corr == FALSE ){
        if( is.matrix(Sigma) && p > 1 ) stop("p need to be equal to the regimes fitted to the tree.")
        if( is.list(Sigma) && !p == length(Sigma) ) stop("length of 'Sigma' need to be equal to the length of 'nu'. So length(Sigma) == length(nu).")
        ## lsig <- length(Sigma)
        lnu <- length(nu)
        if(!p == lnu) stop("length of 'Sigma' and 'nu' need to be equal to the number of R matrices fitted to the data (p). So length(Sigma) == length(nu) == p.")
    }

    ## Maybe add test for the correct parameter values? Would improve user experience.

    ## Note that the user will need to always specify the number of parameters equal to the number of the traits in the model.
    if(den.mu == "unif"){
        mn_trait <- list()
        for( i in 1:r ){
            mn_trait[[i]] <- function(x) dunif(x, min=par.mu[i,1], max=par.mu[i,2], log=TRUE)
        }
    }
    if(den.mu == "norm"){
        mn_trait <- list()
        for( i in 1:r ){
            mn_trait[[i]] <- function(x) dnorm(x, mean=par.mu[i,1], sd=par.mu[i,2], log=TRUE)
        }
    }
    mn <- function(x) sum( sapply(1:r, function(i) mn_trait[[i]](x[i]) ) )

    if( p == 1 ){
        if(den.sd == "unif"){
            sd <- function(x) sum( dunif(x, min=par.sd[1], max=par.sd[2], log=TRUE) )
        }
        if(den.sd == "lnorm"){
            sd <- function(x) sum( dlnorm(x, meanlog=par.sd[1], sdlog=par.sd[2], log=TRUE) )
        }
    } else{

        if(den.sd == "unif"){
            sd_regime <- list()
            for(i in 1:p){
                sd_regime[[i]] <- function(x) sapply(x, function(y) dunif(y, min=par.sd[i,1], max=par.sd[i,2], log=TRUE))
            }
        }
        if(den.sd == "lnorm"){
            sd_regime <- list()
            for(i in 1:p){
                sd_regime[[i]] <- function(x) sapply(x, function(y) dlnorm(y, meanlog=par.sd[i,1], sdlog=par.sd[i,2], log=TRUE))
            }
        }
        sd <- function(x) sum( sapply(1:p, function(i) sd_regime[[i]](x[[i]]) ) )
    }

    if(unif.corr == TRUE){
        ## When the prior on the correlation is set to uniform then all rate matrices will have the same prior.
        corr <- function(x) sum( sapply(x, function(y) log.diwish(y, v=ncol(y)+1, diag(nrow=ncol(y)) ) ) )
    } else{
        if(p == 1){
            if(!is.matrix(Sigma)) stop(" 'Sigma' needs to be a matrix when p == 1. ")
            if(!is.numeric(nu) || nu < r) stop(" 'nu' need to be a numeric value larger than the dimension of 'Sigma'. ")
            center <- (nu - ncol(Sigma) -1) * Sigma
            corr <- function(x) log.diwish(x, v=nu, center)
        } else{
            if(!all(sapply(Sigma, is.matrix))) stop(" 'Sigma' needs to be a list of matrices with length == p. ")
            if(!all(sapply(nu, is.numeric)) || all(nu < r) ) stop(" 'nu' needs to be a vector with numeric values larger than r. ")
            corr_regime <- list()
            center <- list()
            for( i in 1:p){
                center[[i]] <- (nu[i] - ncol(Sigma[[i]]) -1) * Sigma[[i]]
                ## corr_regime[[i]] <- function(x) sapply(x, function(y) log.diwish(y, v=nu[i], center[[i]]) )
                corr_regime[[i]] <- function(x) log.diwish(x, v=nu[i], center[[i]])
            }
            corr <- function(x) sum( sapply(1:p, function(i) corr_regime[[i]](x[[i]]) ) )
        }
    }
    
    return( list(mn, corr, sd) )
}
