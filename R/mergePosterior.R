##' Join two or more independent MCMC chains from the same data and phylogenetic trees by appending them together into a single chain.
##'
##' 
##' @title Merge posterior distributions
##' @param ... any number of posterior distributions as produced by the function 'readMCMC'.
##' @return A merged posterior distribution in the same format.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @examples
##' \donttest{
##' data( centrarchidae )
##' ## Set the limits of the uniform prior on the root based on the observed traits
##' data.range <- t( apply( centrarchidae$data, 2, range ) )
##' ## The step size for the root value can be set given the range we need to sample from:
##' w_mu <- ( data.range[,2] - data.range[,1] ) / 10
##' ## Set a reasonable value for the uniform prior distribution for the standard deviation.
##' ## Here the minimum rate for the traits is 0 and the maximum is 10 ( using 'sqrt(10)' to 
##' ##      transform to standard deviation).
##' par.sd <- cbind(c(0,0), sqrt( c(10,10) ))
##' ## The proposal step on the standard deviation is a multiplier. So 0.2 is good enough 
##' ##       for most cases.
##' w_sd <- matrix(0.2, ncol = 2, nrow = 2)
##' prior <- makePrior(r = 2, p = 2, den.mu = "unif", par.mu = data.range, den.sd = "unif"
##'                    , par.sd = par.sd)
##' ## Run multiple MCMC chains.
##' handle.list <- lapply(1:4, function(x) ratematrixMCMC(data=centrarchidae$data
##'                       , phy=centrarchidae$phy.map, prior=prior, gen=10000
##'                       , w_mu=w_mu, w_sd=w_sd, dir=tempdir()) )
##' ## Read all to a list
##' posterior.list <- lapply(handle.list, readMCMC)
##' ## Merge all posteriors in the list.
##' merged.four <- mergePosterior(posterior.list)
##' ## Merge some of the posteriors.
##' merged.two <- mergePosterior(posterior.list[[1]], posterior.list[[3]])
##' }
mergePosterior <- function(...){
    
    chains <- list(...)

    if( inherits(chains[[1]], what = c("ratematrix_multi_chain", "ratematrix_single_chain") ) ){
        ## Probably a single posterior as input.
        if( length( chains[[1]] ) == 1 ) stop( "Need two or more posterior distributions to merge." )
    } else if( inherits(chains[[1]], what = "list") ){
        ## Input likely a list of posteriors.
        if( inherits(chains[[1]][[1]], what = c("ratematrix_multi_chain", "ratematrix_single_chain") ) ){
            ## Yes, this is a list of posteriors:
            if( length( chains[[1]] ) == 1 ) stop( "Need two or more posterior distributions to merge." )
            class.post <- sapply(chains[[1]], function(x) inherits(x, what = c("ratematrix_multi_chain", "ratematrix_single_chain") ))
            if( !all( class.post ) ) stop( "All input objects need to be posterior distributions." )
            chains <- chains[[1]] ## Decrease one level in the listing structure.
        } else{
            stop( "All input objects need to be posterior distributions." )
        }
    } else{
        ## Don't know what type this is. Returns error.
        stop( "All input objects need to be posterior distributions." )
    }    

    ## Test if all the posteriors provided are of the same class.
    ## Then calculate the value of p, if necessary and proceed.
    ## The type of convergence test will be given by the class of the posterior.
    
    if( sum( sapply(chains, function(x) inherits(x, what=c("ratematrix_single_chain", "ratematrix_multi_chain")) ) ) == length(chains) ){
        
        if( length( unique( sapply(chains, function(x) class(x) ) ) ) > 1 ) stop("Posterior chains need to belong to the same model. \n")
        
        if( length(chains) == 1 ){
            if( inherits(chains[[1]], what=c("ratematrix_single_chain")) ){
                p <- 1
            } else{
                p <- length( chains[[1]]$matrix )
                regimes <- names( chains[[1]]$matrix )
            }
        }
        
        if( length(chains) > 1 ){
            if( inherits(chains[[1]], what=c("ratematrix_single_chain")) ){
                p <- 1
            } else{
                p <- length( chains[[1]]$matrix )
                regimes <- names( chains[[1]]$matrix )
            }
        }
        
    } else{
        stop("Arguments need to be output of 'readMCMC' function. Of class 'ratematrix_single_chain' or 'ratematrix_multi_chain'. \n")
    }
    
    mcmc.join <- list() ## The output of the function.
    ## Need to add element in this order: root, matrix, log.lik.
    mcmc.join$root <- do.call(rbind, lapply(chains, function(x) x$root) )

    if( p == 1 ){
        mcmc.join$matrix <- do.call(c, lapply(chains, function(x) x$matrix) )
        class( mcmc.join ) <- "ratematrix_single_chain"
    } else{        
        mcmc.join$matrix <- lapply(1:p, function(y) do.call(c, lapply(chains, function(x) x$matrix[[y]]) ) )
        names( mcmc.join$matrix ) <- regimes
        class( mcmc.join ) <- "ratematrix_multi_chain"
    }

    return( mcmc.join )
}
