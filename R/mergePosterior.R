##' Join two or more independent MCMC chains from the same data and phylogenetic trees by appending them together into a single chain.
##'
##' 
##' @title Merge posterior distributions
##' @param ... Any number of posterior distributions as produced by the function 'readMCMC'.
##' @return A merged posterior distribution in the same format.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @examples
##' ## Run two (very short) chains and merge results.
##' data(centrarchidae)
##' handle1 <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=1000)
##' posterior1 <- readMCMC(handle1, burn=0.25, thin=1)
##' handle2 <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=1000)
##' posterior2 <- readMCMC(handle2, burn=0.25, thin=1)
##' merged <- mergePosterior(posterior1, posterior2)
##' merged ## Summary of the chain
##' plotRatematrix(merged)
mergePosterior <- function(...){
    
    chains <- list(...)

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
            }
        }
        
        if( length(chains) > 1 ){
            if( inherits(chains[[1]], what=c("ratematrix_single_chain")) ){
                p <- 1
            } else{
                p <- length( chains[[1]]$matrix )
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
        class( mcmc.join ) <- "ratematrix_multi_chain"
    }

    return( mcmc.join )
}
