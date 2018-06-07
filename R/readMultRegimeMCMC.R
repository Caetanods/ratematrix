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
##' @importFrom readr read_lines
##' @noRd
readMultRegimeMCMC <- function(out, burn = 0.5, thin = 1, dir=NULL){

    if(is.null(dir)){
        direct <- out$dir
    } else{
        direct <- dir
    }

    ## In this version the posterior is in a single file.
    mcmc <- read_lines( file=file.path(direct, paste(out$outname, ".", out$ID, ".mcmc", sep="")) )
    header <- mcmc[1]
    header <- as.character( strsplit(x=header, split=";", fixed=TRUE)[[1]] )

    ## Apply thinning and burnin.
    obs.gen <- length( mcmc )-1 ## Compute observed number of samples (Fix for unfinished chains.) Note that header is first line.
    ## Check if we can read the MCMC. In cases when MCMC failed or something.
    if( obs.gen <= thin ) stop("Length of mcmc samples is lower than 'thin'.")
    if( burn <= 0 ){
        post <- seq(2, obs.gen+1, by=thin) ## First line is the header.
    } else{
        post <- seq(round(obs.gen * burn)+1, obs.gen+1, by=thin) ## First line is the header.
    }
    mcmc <- mcmc[post]

    ## Parse the posterior.
    mcmc <- sapply(mcmc, function(x) as.numeric( strsplit(x=x, split=";", fixed=TRUE)[[1]] )
                 , USE.NAMES=FALSE)
    
    ## Define the columns correspondent to the matrix:
    init <- seq(from=1, to=(out$k^2)*out$p, by=out$k^2)
    end <- rev( seq(from=(out$k^2)*out$p, to=1, by=-out$k^2) )
    RR <- list()
    for( i in 1:out$p ){
        RR[[i]] <- lapply(2:dim(mcmc)[2], function(x) matrix( as.matrix(mcmc)[init[i]:end[i],x], nrow=out$k ) )
    }
    ## names(RR) <- c("regime1","regime2") ## This need to be the name of the regimes!!
    
    ## Now find the root values.
    root <- t(sapply(2:dim(mcmc)[2], function(x) as.numeric( as.matrix(mcmc)[(end[out$p]+1):(end[out$p]+out$k),x])))
    colnames(root) <- out$names ## This need to be the names from the data matrix.
    
    ## The the loglik:
    ## lik <- sapply(2:dim(mcmc)[2], function(x) as.matrix(mcmc)[dim(mcmc)[1],x])

    ## out <- list(root = root, matrix = RR, log.lik = lik)
    out <- list(root = root, matrix = RR)
    class(out) <- "ratematrix_multi_chain"
    
    return( out )
}
