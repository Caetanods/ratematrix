##' Read the results of the MCMC for a single evolutionary rate (R) matrix.
##'
##' Function will use 'readr' package to read the output files produced by the Markov chain Monte Carlo
##'    analysis made with the function 'single.R.iwish.mcmc'.
##' @title Read output files of MCMC with a single R matrix.
##' @param out list. The list object returned by the function 'single.R.iwish.mcmc'. See more details in 'help(single.R.iwish.mcmc)'.
##' @param burn numeric. The proportion of the burnin to be pruned from the MCMC chain.
##' @param thin numeric. The thinning of the posterior distribution. Same format as the argument 'by' of the function 'seq'.
##' @param dir string. Directory where output files from the MCMC run are stored. If 'NULL', then function tries to read the files from the current directory. This can be set to a directory different from the 'dir' argument of the 'single.R.iwish.mcmc' function.
##' @return List with the MCMC chain for the phylogenetic mean (root value) and evolutionary rate matrices (R). 'root' are the values for the phylogenetic mean; 'matrix' is a list with the chain of R matrices; 'log.lik' is the log likelihood (not posterior) for the chain.
##' @importFrom readr read_lines
##' @noRd
readSingleRegimeMCMC <- function(out, burn = 0.5, thin = 1, dir=NULL){
    
    if(is.null(dir)){
        direct <- out$dir
    } else{
        direct <- dir
    }

    if( burn == 0 || burn < 0 ){
        post <- seq(round(out$gen)+1, out$gen+1, by=thin)
    } else{
        post <- seq(round(out$gen * burn)+1, out$gen+1, by=thin) ## First line is the header.
    }

    ## In this version the posterior is in a single file.
    mcmc <- read_lines( file=file.path(direct, paste(out$outname, ".", out$ID, ".mcmc", sep="")) )
    header <- mcmc[1]
    header <- as.character( strsplit(x=header, split=";", fixed=TRUE)[[1]] )
    mcmc <- mcmc[post]
    mcmc <- sapply(mcmc, function(x) as.numeric( strsplit(x=x, split=";", fixed=TRUE)[[1]] )
                 , USE.NAMES=FALSE)
    
    ## Define the columns correspondent to the matrix:
    init <- 1
    end <- out$k^2
    RR <- lapply(1:dim(mcmc)[2], function(x) matrix( as.matrix(mcmc)[init:end,x], nrow=out$k ) )
    ## names(RR) <- c("regime1","regime2") ## This need to be the name of the regimes!!
    
    ## Now find the root values.
    root <- t(sapply(1:dim(mcmc)[2], function(x) as.numeric( as.matrix(mcmc)[(end+1):(end+out$k),x])))
    colnames(root) <- out$names ## This need to be the names from the data matrix.
    
    ## The the loglik:
    ## lik <- sapply(1:dim(mcmc)[2], function(x) as.matrix(mcmc)[dim(mcmc)[1],x])

    ## out <- list(root = root, matrix = RR, log.lik = lik)
    out <- list(root = root, matrix = RR)
    class(out) <- "ratematrix_single_chain"
    
    return( out )

}
