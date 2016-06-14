##' Write the MCMC to files.
##'
##' This function writes the chain stored in the 'cache.chain' list of parameters in the MCMC. Then the function will keep only the last generation in the chain and exclude the rest to save space in the RAM memory.
##' @title Write the MCMC to files.
##' @param files list of the files for each of the parameters.
##' @param cache.chain The object with the MCMC generations.
##' @param p the number of matrices fitted in the model.
##' @param chunk The interval of generations when this function should be called to write parameter states to file and clean memory.
##' @return Return a updated version of cache.chain. This chain will now only have the last state in the chain. Function writes to file the content in 'cache.chain'. This function will append the new generations to the previous chunk.
write.to.multfile <- function(files, cache.chain, p, chunk){
    ## Write cache.chain to files.
    ## The number of files will depend on the number of matrices.
    ## files = A list of the open connections.
    ## cache.chain = the current cache.chain
    ## p = the number of R matrices fitted to the data.
    ##ll <- length(cache.chain$chain)
    ll <- chunk+1
    ## Never write the first element. Keep the last element.
    sapply(2:ll, function(x) cat(c(cache.chain$lik[x]),"\n",sep=";",file=files[[1]], append=TRUE))
    sapply(2:ll, function(x) cat(c(cache.chain$chain[[x]][[1]]),"\n",sep=";",file=files[[2]], append=TRUE))
    for(i in 1:p){
        ## Note that 1st position of files is the posteior of the root.
        sapply(2:ll, function(x) cat(c(cache.chain$chain[[x]][[4]][[i]]),"\n",sep=";"
                                    , file=files[[i+2]], append=TRUE))
    }

    ## Keep the last elemente and clear the rest for memory saving.
    keep.chain <- cache.chain$chain[[ll]]
    keep.lik <- cache.chain$lik[ll]
    ##cache.chain$chain <- list()
    cache.chain$chain <- vector(mode="list", length=chunk+1)
    ##cache.chain$lik <- vector()
    cache.chain$lik <- vector(mode="numeric", length=chunk+1)
    cache.chain$chain[[1]] <- keep.chain
    cache.chain$lik[1] <- keep.lik

    return(cache.chain)
}
