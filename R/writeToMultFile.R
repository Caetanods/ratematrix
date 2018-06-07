##' Write the MCMC to files.
##'
##' This function writes the chain stored in the 'cache.chain' list of parameters in the MCMC. Then the function will keep only the last generation in the chain and exclude the rest to save space in the RAM memory.
##' @title Write the MCMC to files.
##' @param files list of the files for each of the parameters.
##' @param cache.chain The object with the MCMC generations.
##' @param p the number of matrices fitted in the model.
##' @param chunk The interval of generations when this function should be called to write parameter states to file and clean memory.
##' @return Return a updated version of cache.chain. This chain will now only have the last state in the chain. Function writes to file the content in 'cache.chain'. This function will append the new generations to the previous chunk.
##' @noRd
writeToMultFile <- function(files, cache.chain, p){
    ## Write cache.chain to files.
    ## The number of files will depend on the number of matrices.
    ## files = A list of the open connections.
    ## cache.chain = the current cache.chain
    ## p = the number of R matrices fitted to the data.
    ## This will write all the output to the same logfile.
    ## First come the p R matrices regimes followed by the root values and then the loglik.

    for( j in 1:p){
        cat( c(cache.chain$chain[[4]][[j]]), sep="; ", file=files[[1]], append=TRUE)
        cat("; ", sep="", file=files[[1]], append=TRUE)
    }
    cat(cache.chain$chain[[1]], sep="; ", file=files[[1]], append=TRUE)
    cat("\n", file=files[[1]], append=TRUE)
}
