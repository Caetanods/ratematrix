#' Extract last sample from a previous MCMC analysis performed using either 'ratematrixMCMC' or 'ratematrixJointMCMC' to use as the starting point state for another MCMC analysis (also using 'ratematrixMCMC' or 'ratematrixJointMCMC', respectivelly). The number of traits and number of regimes need to be the same between the MCMC runs (the phylogeny and the configuration of the regimes can be changed).
#'
#' This function can be used to replicate multiple MCMC runs starting from another MCMC. It can also be used to continue a MCMC chain from its last point, because it was was terminated prematurely or did not converge. Also see 'continueMCMC'. Note that 'continueMCMC' will not work with 'ratematrixJointMCMC' analyses, but this function does.
#' @title Extract sample from MCMC to use as starting point for another MCMC
#' @param handle the output object from the 'ratematrixMCMC' or 'ratematrixJointMCMC' functions.
#' @param dir directory with the output files. If 'NULL' (default), then files are read from the directory chosen when running the MCMC chain using the argument 'dir' of the 'ratematrixMCMC' or 'ratematrixJointMCMC' functions (stored on handle). Otherwise function will read files from 'dir'.
#'
#' @return a sample from the posterior of a previous analysis to be used as the 'start' argument for the 'ratematrixMCMC' or 'ratematrixJointMCMC' functions.
#' @export
#' @author Daniel S. Caetano and Luke J. Harmon
#' @importFrom corpcor decompose.cov
#'
#' @examples
#' \donttest{
#' ## Load data
#' data(anoles)
#' ## Run MCMC. This is just a very short chain.
#' phy <- mergeSimmap(anoles$phy[[1]], drop.regimes = TRUE) ## Turn simmap into phylo.
#' traits <- anoles$data[,c(1,2)] ## The continuous traits
#' ## The predictor data.
#' pred <- setNames(object = as.character(anoles$data$area), nm = rownames(anoles$data))
#' handle <- ratematrixJointMCMC(data_BM = traits, data_Mk = pred, phy = phy
#'                               , gen = 1000, dir = tempdir())
#' ## Load posterior distribution, make plots and check the log.
#' posterior <- readMCMC(handle, burn=0.25, thin=1)
#' plotRatematrix(posterior)
#' plotRootValue(posterior)
#' ## Start another MCMC from the last sample of the previous one.
#' st_pt <- getStartPointFromPosterior(handle = handle)
#' handle_new <- ratematrixJointMCMC(data_BM = traits, data_Mk = pred, phy = phy
#'                                   , start = st_pt, gen=1000, dir=tempdir())
#' post_new <- readMCMC(handle_new, burn=0.25, thin=1)
#' plotRatematrix(post_new)
#' plotRootValue(post_new)
#' }
getStartPointFromPosterior <- function(handle, dir=NULL){
    if(is.null(dir)){
        dir_tmp <- handle$dir
    } else{
        dir_tmp <- dir
    }
    start <- list()

    ## Will generate a starting point object to be used in another MCMC run of the same type.
    if( inherits(x = handle, what = "ratematrix_single_mcmc") ){
        ## Single regime posterior.
        mcmc_sample <- readMCMC(handle = handle, dir = dir_tmp, burn = 0.90, thin = 1)
        last_id <- nrow( mcmc_sample$root )
        start$root <- mcmc_sample$root[last_id,]
        start$matrix <- mcmc_sample$matrix[[last_id]]
        start$sd <- sqrt( decompose.cov(m = start$R)$v )
        class( start ) <- "ratematrix_prior_sample"
        return( start )
    } else if( inherits(x = handle, what = "ratematrix_multi_mcmc") ){
        ## Mult regime, not joint estimate.
        mcmc_sample <- readMCMC(handle = handle, dir = dir_tmp, burn = 0.90, thin = 1)
        last_id <- nrow( mcmc_sample$root )
        start$root <- mcmc_sample$root[last_id,]
        start$matrix <- list()
        for( i in 1:handle$p ) start$matrix[[i]] <- mcmc_sample$matrix[[i]][[last_id]]
        start$sd <- list()
        for( i in 1:handle$p ) start$sd[[i]] <- sqrt(decompose.cov(m = start$matrix[[i]])$v)
        class( start ) <- "ratematrix_prior_sample"
        return( start )
    } else{
        stop("handle must be the output from 'ratematrixMCMC' or 'ratematrixJointMCMC' functions.")
    }
}

