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
##' data( centrarchidae )
##' dt.range <- t( apply( centrarchidae$data, 2, range ) )
##' ## The step size for the root value can be set given the range we need to sample from:
##' w_mu <- ( dt.range[,2] - dt.range[,1] ) / 10
##' par.sd <- cbind(c(0,0), sqrt( c(10,10) ))
##' prior <- makePrior(r=2, p=2, den.mu="unif", par.mu=dt.range, den.sd="unif", par.sd=par.sd)
##' prior.samples <- samplePrior(n = 1000, prior = prior)
##' start.point <- samplePrior(n=1, prior=prior)
##' ## Plot the prior. Red line shows the sample from the prior that will set the 
##' ##      starting point for the MCMC.
##' plotRatematrix(prior.samples, point.matrix = start.point$matrix, point.color = "red"
##'                , point.wd = 2)
##' plotRootValue(prior.samples)
##' handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, prior=prior
##'                          , gen=10000, w_mu=w_mu, dir=tempdir())
##' posterior <- readMCMC(handle, burn = 0.2, thin = 10)
##' ## Again, here the red line shows the starting point of the MCMC.
##' plotRatematrix( posterior, point.matrix = start.point$matrix, point.color = "red"
##'                , point.wd = 2)
##' plotRootValue(posterior)
##' computeESS(mcmc=posterior, p=2)
##' }
computeESS <- function(mcmc, p){
    ## Will need to make some modifications to work with the new class:
    ## Avoid computing quantities for the root value and compute the quantities for the trait samples.

    if( missing(p) ) stop("Need to specify number of regimes as argument 'p'.")
    if( inherits(mcmc, what=c("ratematrix_single_chain", "ratematrix_multi_chain")) ){
        root <- mcmc$root
    }
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

    if( inherits(mcmc, what=c("ratematrix_single_chain", "ratematrix_multi_chain")) ){
        mcmc.root <- coda::mcmc(root)
        ess.root <- coda::effectiveSize(x=mcmc.root)
    }
    
    mcmc.rate <- lapply(rates, coda::mcmc)
    if( k < 3 ){
        mcmc.corr <- lapply(corr_vec, coda::mcmc)
    } else{
        mcmc.corr <- lapply(corr, coda::mcmc)
    }

    ess.rate <- lapply(mcmc.rate, coda::effectiveSize)
    ess.corr <- lapply(mcmc.corr, coda::effectiveSize)

    if( inherits(mcmc, what="ratematrix_poly_chain") ){
        ## Compute the ESS for the trait samples.
        ## Will divide by trait and by category, tips and ancestral values (Meaning that we will pool the trait samples across species)
        range.ess.trait.tip <- range( apply(mcmc$tip_samples, MARGIN = 2, coda::effectiveSize) )
        range.ess.trait.anc <- range( apply(mcmc$anc_samples, MARGIN = 2, coda::effectiveSize) )
        range.ess.trait.tip <- setNames(range.ess.trait.tip, c("min","max"))
        range.ess.trait.anc <- setNames(range.ess.trait.anc, c("min","max"))
    }

    if( inherits(mcmc, what=c("ratematrix_single_chain", "ratematrix_multi_chain")) ){
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
        }
    }
    
    if( inherits(mcmc, what="ratematrix_poly_chain") ){
        if( p == 1 ){
            cat("ESS rates \n")
            print( ess.rate[[1]] )
            cat("\n")
            cat("ESS correlation \n")
            print( ess.corr[[1]] )
            cat("ESS traits sampled at the tips (range) \n")
            print( range.ess.trait.tip )
            cat("ESS traits sampled at the nodes (range) \n")
            print( range.ess.trait.anc )
            cat("\n")
            res <- list( ESS_rates=ess.rate[[1]], ESS_corr=ess.corr[[1]]
                      , ESS_tip_traits=setNames(range.ess.trait.tip, c("min","max"))
                      , ESS_anc_traits=setNames(range.ess.trait.anc, c("min","max")) )
        }
        if( p > 1 ){
            cat("ESS rates \n")
            print( ess.rate )
            cat("\n")
            cat("ESS correlation \n")
            print( ess.corr )
            cat("ESS traits sampled at the tips (range) \n")
            print( range.ess.trait.tip )
            cat("ESS traits sampled at the nodes (range) \n")
            print( range.ess.trait.anc )
            cat("\n")
            res <- list( ESS_rates=ess.rate, ESS_corr=ess.corr
                      , ESS_tip_traits=setNames(range.ess.trait.tip, c("min","max"))
                      , ESS_anc_traits=setNames(range.ess.trait.anc, c("min","max")) )
        }
    }

    return( res )
}
