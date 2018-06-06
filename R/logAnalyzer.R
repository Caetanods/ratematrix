##' Reads the log file produced by the 'ratematrixMCMC' function. Calculates acceptance ratio and shows the trace plot. Check the function 'computeESS' to compute the Effective Sample Size of the posterior distribution.
##'
##' The log shows the acceptance ratio for the parameters of the model and also for each of the phylogenies provided to the 'ratematrixMCMC' function (if more than one was provided as input). Also see function 'ratematrixMCMC' for a brief discussion about acceptance ratio for the parameters in 'Details'.\cr
##' \cr
##' The acceptance ratio is the frequency in which any proposal step for that parameter was accepted by the MCMC sampler. When this frequency is too high, then proposals are accepted too often, which might decrease the efficiency of the sampler to sample from a wide range of the parameter space (the steps are too short). On the other hand, when the acceptance ratio is too low, then the steps of the sampler propose new values that are often outside of the posterior distribution and are systematically rejected by the sampler. Statisticians often suggest that a good acceptance ratio for a MCMC is something close to '0.24'. Our experience is that acceptance ratios between 0.15 and 0.4 will work just fine. Much lower or higher than this might create mixing problems or be too inefficient.\cr
##' \cr
##' If you provided a list of phylogenies to the MCMC chain, then the sampler will randomly sample one of these phylogenies and use it to compute the likelihood of the model at each step of the MCMC. The pool of trees and/or regime configurations provided effectivelly works as a prior distribution. It is important to note that this is not equivalent to a joint estimation of the comparative model of trait evolution and phylogenetic tree, since the moves proposed by the MCMC chain are restricted to the parameters of the phylogenetic comparative model. Some of the phylogenies provided in the pool might be accepted more than others during the MCMC. When this happens, the acceptance ratio for a given tree, or set of trees, will be relativelly lower when compared to the rest. This means that the information presented in such a tree (or trees) is less represented in the posterior distribution than other trees. If this issue happens, we advise users to investigate whether these trees show a different pattern (potentially biologically informative) when compared to the other set of trees. Additionally, one might also repeat the analysis with these trees in separate in order to check whether parameters estimates are divergent.
##' @title Make analysis of the log file of the MCMC chain
##' @param handle the output object from the 'ratematrixMCMC' function.
##' @param burn the proportion of burn-in. A numeric value between 0 and 1.
##' @param thin the number of generations to skip when reading the posterior distributionfrom the files. Since the files contain each step of the sampler, one can check the posterior with different 'thin' values without the need of reanalyses.
##' @param show.plots whether to show a trace plot of the log-likelihood and the acceptance ratio. Default is TRUE.
##' @param print.result whether to print the results of the acceptance ratio to the screen. Default is TRUE.
##' @param dir the directory where to find the log file. If set to 'NULL' (default), the function will search for the files in the same directory that the MCMC chain was made (stored in handle$dir).
##' @return A named vector with the acceptance ratio for the whole MCMC and each of the parameters of the model. If a list of phylogenetic trees was provided to the MCMC chain, then the output is a list with the acceptance ratio for the parameters and a table showing the frequency in which each of the phylogenies was accepted in a move step.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @examples
##' \donttest{
##' ## Load data
##' data(centrarchidae)
##' ## Run MCMC. This is just a very short chain.
##' handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=10000
##'                          , dir=tempdir())
##' ## Load posterior distribution, make plots and check the log.
##' posterior <- readMCMC(handle, burn=0.1, thin=10)
##' plotRatematrix(posterior)
##' logAnalyzer(handle, burn=0.1, thin=10)
##' }
logAnalyzer <- function(handle, burn=0.25, thin=100, show.plots=TRUE, print.result=TRUE, dir=NULL){
    ## Need to read the log output. Then make some plots or such.

    ## Check the directory.
    if(is.null(dir)){
        direct <- handle$dir
    } else{
        direct <- dir
    }

    ## Read the log and create a matrix.
    log.mcmc <- read_lines( file=file.path(direct, paste(handle$outname, ".", handle$ID, ".log", sep="")) )
    header <- log.mcmc[1]
    header <- as.character( strsplit(x=header, split=";", fixed=TRUE)[[1]] )

    ## Compute the posterior based on the observed samples.
    obs.gen <- length( log.mcmc )-1 ## Compute observed number of samples (Fix for unfinished chains.) Note that header is first line.
    if( burn <= 0 ){
        post <- seq(2, obs.gen+1, by=thin) ## First line is the header.
    } else{
        post <- seq(round(obs.gen * burn)+1, obs.gen+1, by=thin) ## First line is the header.
    }
    log.mcmc <- log.mcmc[post]

    ## Parse the samples.
    log.mcmc <- t( sapply(log.mcmc, function(x) as.numeric( strsplit(x=x, split=";", fixed=TRUE)[[1]] )
                        , USE.NAMES=FALSE) )
    colnames( log.mcmc ) <- header

    ## Make analyses.

    ## The accept ratio for all parameters:
    accept.all <- sum( as.logical( log.mcmc[,1]) ) / nrow(log.mcmc)
    ## Compute the cumulative acceptance ratio for all the parameters across the generations. This is based on the log after burn-in and thinning.
    cum.accept <- cumsum(log.mcmc[,1]) / 1:length(log.mcmc[,1])
    ## The accept ratio for each of the parameters in separate.
    accept.par <- sapply(2:4, function(x) sum( log.mcmc[ as.logical( log.mcmc[,x] ), 1] ) / length( log.mcmc[ as.logical( log.mcmc[,x] ), 1] ) )
    ## Make report object:
    accept <- c(accept.all, accept.par)
    names(accept) <- c("all", "correlation", "sd", "root")
    
    ## The proportion of times that each phylo of the list was present when a proposal was accepted.
    ## This is the 'mixing' for the pool of phylogenetic trees.
    mix.phylo <- table( log.mcmc[ as.logical(log.mcmc[,1]), 5] ) / length( log.mcmc[ as.logical(log.mcmc[,1]), 5 ] )

    ## Create object to help making the plot.
    at.gen <- round( seq(from=1, to=length(post), length.out = 5) )
    labels.gen <- round( seq(from=post[1], to=post[length(post)], length.out = 5) )

    if( print.result==TRUE ){
        cat("Acceptance ratio for the MCMC and parameters:\n")
        print( accept )
        if( is.list(handle$phy[[1]]) ){
            cat("\n")
            cat("Acceptance ratio for each phylogeny:\n")
            print( mix.phylo )
        }
    }

    if( show.plots==TRUE ){
        ## Log-likelihood and acceptance ration trace plot:
        old.par <- graphics::par(no.readonly = TRUE)

        graphics::par( mfrow = c(2,1) )
        graphics::par(mar = c(1, 0, 0, 0), oma = c(3, 4, 2, 0))
        graphics::plot(x=1:nrow( log.mcmc ), y=log.mcmc[,6], type="l", axes=FALSE, xlab="", ylab="")
        graphics::axis(side=2)
        graphics::mtext("Log-likelihood", side=2, line=2, cex=1)
        graphics::plot(x=1:nrow( log.mcmc ), y=cum.accept, type="l", axes=FALSE, xlab="", ylab="")
        graphics::axis(side=2)
        graphics::axis(side=1, at=at.gen, labels=labels.gen)
        graphics::mtext("Generations", side=1, line=2, cex=1)
        graphics::mtext("Acceptance ratio", side=2, line=2, cex=1)

        graphics::par(old.par)
    }

    if( is.list(handle$phy[[1]]) ){
        return( list(accept.ratio=accept, mix.phylo=mix.phylo) )
    } else{
        return( accept.ratio=accept )
    }
}
