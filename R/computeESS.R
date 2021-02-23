##' Computes the Effective Sample Size (ESS) for the parameters of the model from the MCMC samples.
##'
##' Function uses 'coda' function 'effectiveSize' to compute the ESS for each of the parameters of the model separatelly. Values for the ESS is too low indicates poor mixing for the parameter of the model.
##' @title Compute the ESS for the MCMC samples
##' @param mcmc Posterior distribution object. Same as output from 'readMCMC' function.
##' @param p Number of evolutionary rate matrix regimes fitted to the phylogenetic tree.
##' @return A list object with the ESS value for the root, evolutionary rates, and evolutionary correlations among the traits.
##' @author Daniel Caetano and Luke Harmon
##' @export
##' @importFrom coda mcmc
##' @importFrom coda effectiveSize
##' @examples
##' \donttest{
##' library( ratematrix )
##' data( centrarchidae )
##' ## Run multiple MCMC chains.
##' handle.list <- lapply(1:4, function(x) ratematrixMCMC(data=centrarchidae$data
##'                       , phy=centrarchidae$phy.map, gen=10000, dir=tempdir()) )
##' ## Read all to a list
##' posterior.list <- lapply(handle.list, function(x) readMCMC(handle = x, burn = 0.5))
##' ## Merge all posteriors.
##' merged.four <- mergePosterior(posterior.list)
##' ## Compute the ESS for the merged posterior:
##' ess <- computeESS(mcmc=merged.four, p=2)
##' }
computeESS <- function(mcmc, p){
    if( missing(p) ) stop("Need to specify number of regimes as argument 'p'.")
    root <- mcmc$root
    rates <- list()
    corr <- list()

    if( p == 1 ){

        ## Check the accuracy of the p argument.
        if( !is.matrix( mcmc$matrix[[1]] ) ) stop( "Number of regimes might be wrong." )

        k <- ncol(mcmc$matrix[[1]])
        
        rates[[1]] <- t( sapply(mcmc$matrix, function(x) diag(x) ) )
        upper <- upper.tri(mcmc$matrix[[1]])
        corr[[1]] <- t( sapply(mcmc$matrix, function(x) c( stats::cov2cor(x)[upper] ) ) )
        if( k < 3 ){
            corr_vec <- list()
            corr_vec[[1]] <- corr[[1]][1,]
        }
    }

    if( p > 1){

        ## Check the accuracy of the p argument.
        if( !p == length( mcmc$matrix ) ) stop( "Number of regimes might be wrong." )
        check.p <- sapply(1:p, function(x) is.matrix( mcmc$matrix[[x]][[1]] ) )
        if( !all(check.p) ) stop( "Number of regimes might be wrong." )
        
        k <- ncol( mcmc$matrix[[1]][[1]] )
        
        for( i in 1:p ){
            rates[[i]] <- t( sapply(mcmc$matrix[[i]], function(x) diag(x) ) )
            upper <- upper.tri(mcmc$matrix[[i]][[1]])
            corr[[i]] <- t( sapply(mcmc$matrix[[i]], function(x) c( stats::cov2cor(x)[upper] ) ) )
        }
        if( k < 3 ){
            ## Elements of 'corr' will be vectors.
            corr_vec <- list()
            corr_vec <- lapply(1:p, function(x) corr[[x]][1,])
        }
    }

    mcmc.root <- coda::mcmc(root)
    mcmc.rate <- lapply(rates, coda::mcmc)
    if( k < 3 ){
        mcmc.corr <- lapply(corr_vec, coda::mcmc)
    } else{
        mcmc.corr <- lapply(corr, coda::mcmc)
    }

    ess.root <- coda::effectiveSize(x=mcmc.root)
    ess.rate <- lapply(mcmc.rate, coda::effectiveSize)
    ess.corr <- lapply(mcmc.corr, coda::effectiveSize)
    
    if( p == 1 ){
        cat("ESS root values \n")
        print( ess.root )
        cat("\n")
        cat("ESS rates \n")
        print( ess.rate[[1]] )
        cat("\n")
        cat("ESS correlation \n")
        print( ess.corr[[1]] )
        cat("\n")
        res <- list( ess.root, ess.rate[[1]], ess.corr[[1]] )
        names( res ) <- c("ESS_root","ESS_rates","ESS_corr")
        return( res )
    }

    if( p > 1 ){
        cat("ESS root values \n")
        print( ess.root )
        cat("\n")
        cat("ESS rates \n")
        print( ess.rate )
        cat("\n")
        cat("ESS correlation \n")
        print( ess.corr )
        cat("\n")
        res <- list( ess.root, ess.rate, ess.corr )
        names( res ) <- c("ESS_root","ESS_rates","ESS_corr")
        return( res )
    }
    
}
