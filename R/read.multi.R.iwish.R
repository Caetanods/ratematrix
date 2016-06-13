##' Read the results of the MCMC for a two or more evolutionary rate (R) matrices.
##'
##' Function will use 'readr' package to read the output files produced by the Markov chain Monte Carlo
##'    analysis made with the function 'multi.R.iwish.mcmc'.
##' @title Read output files of MCMC with multiple R matrices.
##' @param out list. The list object returned by the function 'multi.R.iwish.mcmc'. See more details in 'help(multi.R.iwish.mcmc)'.
##' @param burn numeric. The proportion of the burnin to be pruned from the MCMC chain.
##' @param thin numeric. The thinning of the posterior distribution. Same format as the argument 'by' of the function 'seq'.
##' @param dir string. Directory where output files from the MCMC run are stored. If 'NULL', then function tries to read the files from the current directory. This can be set to a directory different from the 'dir' argument of the 'multi.R.iwish.mcmc' function.
##' @return List with the MCMC chain for the phylogenetic mean (root value) and evolutionary rate matrices (R). 'root' are the values for the phylogenetic mean; 'matrix' is a list of length equal to the number of matrices fitted to the tree, each of those are lists with the chain of respective R matrices; 'log.lik' is the log likelihood (not posterior) for the chain.
##' @export
##' @importFrom readr read_lines
read.multi.R.iwish <- function(out, burn = 0.5, thin = 1, dir=NULL){
    ## Function to read the files written by iwish fast.
    ## Depends on the package 'readr'.
    ## Function takes as arguments the output of 'singleR_iwish_fast'.
    ## Optionally a different directory from the original can be provided.
    ## Returns a list of root values and R matrices. Read for other functions.
    ## In the future this function will read the whole output of the MCMC.
    ## out = Output object from 'singleR_iwish_fast'.
    ## burn = The percent of burn-in. Value between 0 and 1. The burn in will be taken out
    ##       while the function read the output. This will save memory and speed up the process
    ##       whitout throwing the data away.
    ## thin = The thinning of the posterior. Same format as the argument 'by' in function 'seq'.
    ## dir = Optional directory to find results. Different from out$dir.

    if(is.null(dir)){
        direct <- out$dir
    } else{
        direct <- dir
    }

    post <- seq(round(out$gen * burn), out$gen, by=thin)

    lik <- read_lines( file=file.path(direct, paste(out$outname, ".", out$ID, ".loglik", sep="")) )
    lik <- lik[post]
    lik <- sapply(lik, function(x) as.numeric( strsplit(x=x, split=";", fixed=TRUE)[[1]] )
                  , USE.NAMES=FALSE)
    
    MM <- list() ## List object to keep the list of matrices.
    for(i in 1:out$p){
        RR <- read_lines( file=file.path(direct, paste(out$outname,".",out$ID,".",i,".matrix", sep="")) )
        RR <- RR[post]
        RR <- t( sapply(RR, function(x) as.numeric( strsplit(x=x, split=";", fixed=TRUE)[[1]] )
                      , USE.NAMES=FALSE) )
        MM[[i]] <- lapply(1:dim(RR)[1], function(x) (matrix( as.matrix(RR)[x,], nrow=out$k) ) )
    }
    rm(RR) ## To save memory.
    
    root <- read_lines( file=file.path(direct, paste(out$outname,".", out$ID, ".root", sep="")) )
    root <- root[post]
    root <- t( sapply(root, function(x) as.numeric( strsplit(x=x, split=";", fixed=TRUE)[[1]] )
                  , USE.NAMES=FALSE) )
    colnames(root) <- out$names
        
    return( list(root = root, matrix = MM, log.lik = lik) )
}
