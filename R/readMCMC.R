##' Reads the output files from the MCMC with the posterior distribution of the chains.
##' 
##' @title Read the MCMC output files
##' @param handle the output object from the 'ratematrixMCMC' function.
##' @param burn the proportion of the burnin to be pruned from the MCMC chain. A number between 0 and 1 (default is 0).
##' @param thin the thinning of the posterior distribution. A number with the interval between each MCMC step to be kept in the posterior distribution (default is 1).
##' @param dir directory with the output files. If 'NULL' (default), then files are read from the directory chosen when running the MCMC chain using the argument 'dir' of the 'ratematrixMCMC' function (stored on handle). Otherwise function will read files from 'dir'.
##' @return List with the MCMC chain for the phylogenetic mean (root value) and evolutionary rate matrices (R). *root* are the values for the phylogenetic mean in matrix format; *matrix* is a list of length equal to the number of rate regimes fitted to the tree, each of those are lists with the chain of respective R matrices.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @examples
##' \donttest{
##' ## Load data
##' data(centrarchidae)
##' ## Run MCMC. This is just a very short chain.
##' handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=1000
##'                          , dir=tempdir())
##' ## Load posterior distribution, make plots and check the log.
##' posterior <- readMCMC(handle, burn=0.25, thin=1)
##' plotRatematrix(posterior)
##' plotRootValue(posterior)
##' plotPrior(handle)
##' plotPrior(handle, root=TRUE)
##' logAnalyzer(handle)
##' }
readMCMC <- function(handle, burn=0, thin=1, dir=NULL){
    
    if( !inherits(handle, what=c("ratematrix_single_mcmc", "ratematrix_multi_mcmc")) ){
        stop( "Argument 'handle' need to be the output of the 'ratematrixMCMC' function." )
    }
    
    if(handle$p == 1){
        mcmc <- readSingleRegimeMCMC(out=handle, burn=burn, thin=thin, dir=dir)
        colnames( mcmc$root ) <- handle$trait.names
    }
    
    if(handle$p > 1){
        mcmc <- readMultRegimeMCMC(out=handle, burn=burn, thin=thin, dir=dir)
        colnames( mcmc$root ) <- handle$trait.names
        names( mcmc$matrix ) <- handle$regime.names
    }
    
    return( mcmc )
}
