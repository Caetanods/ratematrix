##' Merge two or more posterior distributions
##'
##' Merge two or more posterior distributions based in the same dataset. If you have multiple MCMC chains for the same dataset you can merge them.
##' @title Merge posteriors
##' @param ... Any number of posterior distributions as produced by the function 'read.multi.R.iwish' or 'read.single.R.iwish'
##' @param p Number of R matrices (regimes) fitted to the tree.
##' @return A merged posterior distribution in the same format.
##' @export
mergePosterior <- function(..., p){
    ## Will take mcmc outputs as the arguments. Will return a MCMC output with the merge of all imputs.
    mcmc.list <- list(...) ## This will make a list with all the imput arguments.
    
    mcmc.join <- list() ## The output of the function.
    ## Need to add element in this order: root, matrix, log.lik.
    mcmc.join$root <- do.call(rbind, lapply(mcmc.list, function(x) x$root) )

    if( p == 1 ){
        mcmc.join$matrix <- do.call(c, lapply(mcmc.list, function(x) x$matrix) )
    } else{        
        mcmc.join$matrix <- lapply(1:p, function(y) do.call(c, lapply(mcmc.list, function(x) x$matrix[[y]]) ) )
    }

    mcmc.join$log.lik <- do.call(c, lapply(mcmc.list, function(x) x$log.lik ) )
    class( mcmc.join ) <- "ratematrix_merged_posterior"
    return( mcmc.join )
}
