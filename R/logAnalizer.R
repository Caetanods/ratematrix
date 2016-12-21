##' Reads the log file produced by the 'ratematrixMCMC' function. Calculates acceptance ratio and shows the trace plot.
##'
##' The log shows the acceptance ratio for the parameters of the model and also for each of the phylogenies provided to the 'ratematrixMCMC' function, if more than one was provided as imput. The function 'ratematrixMCMC' also provides a brief discussion about acceptance ratio for the parameters in the 'Details'.\cr
##' \cr
##' The acceptance ratio is the frequency in which any proposal step for that parameter was accepted by the MCMC sampler. When this frequency is too high, then proposals are accepted too often, which might decrease the efficiency of the sampler to sample from a wide range of the parameter space (the steps are too short). On the other hand, when the acceptance ratio is too low, then the steps of the sampler propose new values that are often outside the posterior distribution and are systematically rejected by the sampler. Statisticians often suggest that a good acceptance ratio for a MCMC is something close to '0.24'. \cr
##' \cr
##' If you provided a list of phylogenies to the MCMC chain, then the sampler will randomly sample one of these phylogenies when evaluating the likelihood of the model. Some of those phylogenies might provide very distinct values of log-likelihood for the model (given the same vector of parameter values). Thus, it is possible that an 'unlikely' phylogeny given the data and the model or an 'unlikely' regime configuration for the evolutionary rate matrices might 'get stuck' during the MCMC. If this happen, then everytime the sampler chooses this phylogeny, the proposal will be rejected regardless of the parameter values for the model. This issue does not affect the convergence of the parameter values for the model, however the posterior distribution will not incorporate the phylogenetic information provided by this particular phylogeny. One way to spot if any of the phylogenies is showing this problem is to look to the acceptance ratio for each of the phylogenies in the table provided by this function. If the MCMC is integrating well every phylogeny provided in this list, then we would expect similar acceptance ratios acrooss the list of phylogenies. Those with significative lower acceptance ratio might represent very unlikely phylogenies and/or rate regime configurations. There is no requirement for this acceptance ratio to be ~ 0.24, just that the value need to be similar across the different phylogenies in the list.
##' @title Make analysis of the log file of the MCMC chain
##' @param handle the output object from the 'ratematrixMCMC' function.
##' @param burn the proportion of burn-in. A numeric value between 0 and 1.
##' @param thin the number of generations to skip when reading the posterior distribution.
##' @param show.plots whether to show a trace plot of the log-likelihood and the acceptance ratio.
##' @param print.result whether to print the results of the acceptance ratio to the screen.
##' @param dir the directory where to find the log file. If set to 'NULL' (default), the function will search for the files in the same directory that the MCMC chain was made (handle$dir).
##' @return A named vector with the acceptance ratio for the whole MCMC and each of the parameters of the model. If a list of phylogenetic trees was provided to the MCMC chain, then the output is a list with the acceptance ratio for the parameters and a table showing the frequency in which each of the phylogenies was accepted in a move step.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @seealso \code{\link{ estimateTimeMCMC }} to estimate the time for the MCMC chain, \code{\link{ readMCMC }} for reading the output files, \code{\link{ plotPrior }} for plotting the prior, \code{\link{ plotRatematrix }} and \code{\link{ plotRootValue }} for plotting the posterior,  \code{\link{ checkConvergence }} to check convergence, \code{\link{ testRatematrix }} to perform tests, and \code{\link{ logAnalizer }} to read and analyze the log file.
##' @examples
##' \donttest{
##' ## Load data
##' data(centrarchidae)
##' ## Run MCMC. This is just a very short chain.
##' handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=1000)
##' ## Load posterior distribution, make plots and check the log.
##' posterior <- readMCMC(handle, burn=0.25, thin=1)
##' plotRatematrix(posterior)
##' plotRootValue(posterior)
##' plotPrior(handle)
##' logAnalizer(handle, burn=0.25, thin=1)
##' }
logAnalizer <- function(handle, burn=0.25, thin=100, show.plots=TRUE, print.result=TRUE, dir=NULL){
    ## Need to read the log output. Then make some plots or such.

    ## Check the directory.
    if(is.null(dir)){
        direct <- handle$dir
    } else{
        direct <- dir
    }

    ## Read the log and create a matrix.
    post <- seq(round(handle$gen * burn)+1, handle$gen+1, by=thin) ## First line is the header.
    log.mcmc <- read_lines( file=file.path(direct, paste(handle$outname, ".", handle$ID, ".log", sep="")) )
    header <- log.mcmc[1]
    header <- as.character( strsplit(x=header, split=";", fixed=TRUE)[[1]] )
    log.mcmc <- log.mcmc[post]
    log.mcmc <- t( sapply(log.mcmc, function(x) as.numeric( strsplit(x=x, split=";", fixed=TRUE)[[1]] )
                        , USE.NAMES=FALSE) )
    colnames( log.mcmc ) <- header

    ## Make analyses:
    accept <- sapply(1:4, function(x) sum(log.mcmc[,x]) / dim(log.mcmc)[1] )
    names(accept) <- c("All", "correlation", "sd", "root")

    mix.phylo <- table( log.mcmc[ which(log.mcmc[,1] == 1), 5] ) / length( log.mcmc[ log.mcmc[,1],5 ] )

    at.gen <- round( seq(from=1, to=length(post), length.out = 5) )
    labels.gen <- round( seq(from=post[1], to=post[length(post)], length.out = 5) )
    cum.accept <- cumsum(log.mcmc[,1]) / 1:length(log.mcmc[,1])

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
        old.par <- par(no.readonly = TRUE)

        par( mfrow = c(2,1) )
        par(mar = c(1, 0, 0, 0), oma = c(3, 4, 2, 0))
        plot(x=1:dim( log.mcmc )[1], y=log.mcmc[,6], type="l", axes=FALSE, xlab="", ylab="")
        axis(side=2)
        mtext("Log-likelihood", side=2, line=2, cex=1)
        plot(x=1:dim( log.mcmc )[1], y=cum.accept, type="l", axes=FALSE, xlab="", ylab="")
        axis(side=2)
        axis(side=1, at=at.gen, labels=labels.gen)
        mtext("Generations", side=1, line=2, cex=1)
        mtext("Acceptance ratio", side=2, line=2, cex=1)

        par(old.par)
    }

    if( is.list(handle$phy[[1]]) ){
        return( list(accept.ratio=accept, mix.phylo=mix.phylo) )
    } else{
        return( accept.ratio=accept )
    }
}
