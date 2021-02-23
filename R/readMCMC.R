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
##' handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=10000
##'                          , dir=tempdir())
##' ## Load posterior distribution.
##' ## We can set the burn-in value and the thinning.
##' posterior <- readMCMC(handle, burn=0.25, thin=1)
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

##' Read the results of the MCMC chain for the predictor regimes.
##'
##' Function will use 'readr' package to read the output files produced by the Markov chain Monte Carlo
##'    analysis focused on the transition rates between the predictor states (the Mk model) made as part
##'    of the 'ratematrixJointMCMC' analysis.
##' @title Read output files of MCMC for the predictor regimes.
##'
##' @param handle the output of the 'ratematrixJointMCMC' function.
##' @param burn the proportion of the burnin to be pruned from the MCMC chain. A number between 0 and 1 (default is 0).
##' @param thin the thinning of the posterior distribution. A number with the interval between each MCMC step to be kept in the posterior distribution (default is 1).
##' @param dir directory with the output files. If 'NULL' (default), then files are read from the directory chosen when running the MCMC chain using the argument 'dir' of the 'ratematrixJointMCMC' function (stored on handle). Otherwise function will read files from 'dir'.
##' @param return_Simmap if the output should be stochastic mapping simulations conditioned on a sample of Q matrices from the posterior distribution (TRUE) or a list of sampled Q matrices.
##' @param nsims number of stochastic mapping simulations to be done. These are based on the last 'nsims' sampled by the MCMC.
##' @param max_nshifts the maximum number of state transitions in any branch of the phylogeny when making stochastic mapping simulations. If you get errors in the stochastic mapping step, try increasing this value. However, this indicates a VERY fast transition rate. See 'fastSimmap' function for more information.
##'
##' @return List with the MCMC chain for the transition matrices or a list of stochastic mappings. See parameter 'return_Simmap'.
##' @importFrom readr read_lines
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
readMCMC_Mk <- function(handle, burn=0, thin=1, dir=NULL, return_Simmap=FALSE, nsims=10, max_nshifts=10000){
    ## Read the Q matrix behind the joint MCMC estimates.
    if(is.null(dir)){
        direct <- handle$dir
    } else{
        direct <- dir
    }

    ## In this version the posterior is in a single file.
    mcmc <- read_lines(file = file.path(direct, paste(handle$outname, handle$ID, "Q_mcmc", sep = ".")) )
    header <- mcmc[1]
    header <- as.character( strsplit(x=header, split=";", fixed=TRUE)[[1]] )

    ## Apply thinning and burn-in.
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

    ## Rebuild the Q matrices.
    get_named_Q <- function(x){
        ## Temporary helper function.
        Q_tmp <- buildQ(vec_Q = x, size = handle$p, model_Q = handle$model_Mk)
        rownames( Q_tmp ) <- handle$regime.names
        colnames( Q_tmp ) <- handle$regime.names
        return( Q_tmp )
    }
    Q_list <- lapply(1:ncol(mcmc), function(w) get_named_Q(x = mcmc[,w]) )

    ## Check if to estimate stochastic maps before returning.
    if( return_Simmap ){
        st_vec <- sapply(1:ncol(handle$data_Mk), function(x) rownames(handle$data_Mk)[as.logical(handle$data_Mk[,x])])
        mk_dt <- setNames(object = st_vec, nm = colnames(handle$data_Mk))
        ## Take the last nsims:
        if( nsims >= length(Q_list) ) stop("'nsims' is larger than the number of posterior samples.")
        Q_simmap <- Q_list[(length(Q_list) - nsims + 1):length(Q_list)]

        ## Check how many phylogenies we have, if we have multiple, then take a sample.
        multi_phylo <- check_phy_list(handle$phy)
        if( multi_phylo ){ ## More than one phylo
            repo <- length(handle$phy) < nsims
            id_phylo <- sample(x = 1:length(handle$phy), size = nsims, replace = repo)
            phy_sample <- handle$phy[id_phylo]
            phy_sample <- lapply(phy_sample, function(x) mergeSimmap(phy = x, drop.regimes = TRUE))
            print(paste0("Making ", nsims, " stochastic mapping simulations."))
            print("Decrease the number of 'nsims' if taking too long.")
            simmap <- lapply(1:nsims, function(x) fastSimmap(tree = phy_sample[[x]], x = mk_dt, Q = Q_simmap[[x]], pi = handle$root_Mk, max_nshifts = 10000))
            return( simmap )
        } else{ ## A single phylogeny
            phy_sample <- mergeSimmap(phy = handle$phy, drop.regimes = TRUE)
            print(paste0("Making ", nsims, " stochastic mapping simulations."))
            print("Decrease the number of 'nsims' if taking too long.")
            simmap <- lapply(1:nsims, function(x) fastSimmap(tree = phy_sample, x = mk_dt, Q = Q_simmap[[x]], pi = handle$root_Mk, max_nshifts = 10000))
            return( simmap )
        }
    } else{
        ## This will create a new class that we can use to print or do other things in the future.
        class(Q_list) <- append(x = class(Q_list), values = "ratematrix_Mk_chain")
        return(Q_list)
    }
}

