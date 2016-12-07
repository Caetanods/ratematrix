##' Read the results of a MCMC chain.
##'
##' Function will use 'readr' package to read the output files produced by the Markov chain Monte Carlo
##'    analysis performed with the function 'ratematrixMCMC'.
##' @title Read the output files from the MCMC.
##' @param out the output object from the 'ratematrixMCMC' function.
##' @param burn the proportion of the burnin to be pruned from the MCMC chain. A number between 0 and 1 (default is 0.25).
##' @param thin the thinning of the posterior distribution. A number with the interval between each MCMC step to be kept in the posterior distribution (default is 100).
##' @param dir directory where output files from the MCMC run are stored. If 'NULL' (default), then files are read from the directory chose when running the MCMC chain using the argument 'dir' of the 'ratematrixMCMC' function. Otherwise function will read files from 'dir'.
##' @return List with the MCMC chain for the phylogenetic mean (root value) and evolutionary rate matrices (R). 'root' are the values for the phylogenetic mean; 'matrix' is a list of length equal to the number of matrices fitted to the tree, each of those are lists with the chain of respective R matrices; 'log.lik' is the log likelihood (not posterior) for the chain.
##' @export
##' @author daniel
readMCMC <- function(out, burn=0.25, thin=100, dir=NULL){
    ## Need to create a way to load the output files without a "out" file. This might be important in the case that the user forgot to save the output file of the MCMC run.
    
    if( !inherits(out, what=c("ratematrix_single_mcmc", "ratematrix_multi_mcmc")) ){
        stop( "Argument 'out' need to be the output of the 'ratematrixMCMC' function." )
    }
    
    if(out$p == 1){
        mcmc <- readSingleRegimeMCMC(out=out, burn=burn, thin=thin, dir=dir)
    }
    if(out$p > 1){
        mcmc <- readMultRegimeMCMC(out=out, burn=burn, thin=thin, dir=dir)
    }
    
    return( mcmc )
}
